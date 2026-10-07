#!/usr/bin/env python3
"""Train the official PRSformer architecture on aligned UKB data.

PRSformer was developed by 23andMe, Inc. This independently written training
adapter imports the user's checkout; it does not redistribute the upstream code.
See https://github.com/23andMe/PRSformer/blob/main/LICENSE.txt.

Adaptation: three masked tasks, continuous-trait residual MSE and binary T2DM
logistic-offset BCE. Covariate baselines and all transforms use training data.
The official genotype encoder and phenotype head remain unchanged.
"""

from __future__ import annotations

import argparse
import ast
import contextlib
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import time
import types

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import rankdata

sys.dont_write_bytecode = True
SUPPORTED_TRAITS = {"height": "continuous", "ldl": "continuous", "t2dm": "binary"}


# 🚩 Input identity, alignment and training-only transforms
def csv_list(value):
	return [item.strip() for item in (value or "").split(",") if item.strip()]


def digest_json(value):
	return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def file_sha256(path):
	hash_value = hashlib.sha256()
	with open(path, "rb") as handle:
		for block in iter(lambda: handle.read(8 * 1024 ** 2), b""):
			hash_value.update(block)
	return hash_value.hexdigest()


def validate_genotypes(values):
	if not np.isfinite(values).all() or not np.all((values == -1) | ((values >= 0) & (values <= 2))):
		raise ValueError("Genotypes must contain dosage 0..2 or the missing value -1; do not standardize SNP dosages.")


def read_inputs(args, traits, require_training=True):
	input_paths = [Path(args.data), Path(args.variants), Path(args.genotypes)]
	initial_stamps = [(path.stat().st_size, path.stat().st_mtime_ns, path.stat().st_ino) for path in input_paths]
	data = pd.read_csv(args.data, sep="\t", dtype={"eid": str, "target": str, "split": str})
	required = ["eid", "target", "split"] + traits
	missing = [column for column in required if column not in data]
	if missing:
		raise ValueError(f"Missing columns in --data: {missing}")
	if data["eid"].isna().any() or data["eid"].duplicated().any():
		raise ValueError("--data must have one unique nonmissing eid per genotype row.")
	if data["target"].isna().any() or data["split"].isna().any():
		raise ValueError("Every sample needs a target label and an explicit split.")
	allowed_splits = {"train", "validation", "test"}
	if not set(data["split"]).issubset(allowed_splits):
		raise ValueError("split must be exactly train, validation, or test.")
	if require_training and set(data["split"]) != allowed_splits:
		raise ValueError("Training requires nonempty train, validation, and test partitions.")
	if not np.any(data["split"] == "test"):
		raise ValueError("No held-out test samples were provided.")
	genotypes = np.load(args.genotypes, mmap_mode="r", allow_pickle=False)
	if genotypes.ndim != 2 or genotypes.shape[0] != len(data):
		raise ValueError("The sample-major genotype matrix must have exactly one row per row of --data, in the same order.")
	if not genotypes.flags.c_contiguous:
		raise ValueError("--genotypes must be a C-contiguous sample-major .npy matrix.")
	if genotypes.dtype not in (np.dtype("int8"), np.dtype("float16"), np.dtype("float32")):
		raise ValueError("Supported genotype dtypes are int8, float16 and float32.")
	variants = pd.read_csv(args.variants, sep="\t", dtype=str)
	variant_columns = ["CHR", "BP", "SNP", "REF", "ALT"]
	if any(column not in variants for column in variant_columns):
		raise ValueError("--variants needs CHR, BP, SNP, REF, ALT columns in genotype-column order.")
	if len(variants) != genotypes.shape[1] or not len(variants):
		raise ValueError("Variant rows must exactly match the genotype columns.")
	if variants[variant_columns].isna().any().any() or variants["SNP"].duplicated().any():
		raise ValueError("Variant metadata must be complete with unique SNP IDs.")
	chromosomes = pd.to_numeric(variants["CHR"].str.replace(r"^chr", "", regex=True, case=False), errors="raise").to_numpy(float)
	positions = pd.to_numeric(variants["BP"], errors="raise").to_numpy(float)
	if not np.all(np.isfinite(chromosomes) & (chromosomes == np.floor(chromosomes)) & (chromosomes >= 1) & (chromosomes <= 22)):
		raise ValueError("PRSformer inputs must use autosomes 1..22.")
	if not np.all(np.isfinite(positions) & (positions == np.floor(positions)) & (positions > 0)):
		raise ValueError("Variant BP must be a positive integer.")
	ordered = (chromosomes[1:] > chromosomes[:-1]) | ((chromosomes[1:] == chromosomes[:-1]) & (positions[1:] >= positions[:-1]))
	if not ordered.all():
		raise ValueError("Input columns must be ordered by numeric chromosome and then base-pair position.")
	if not variants["REF"].str.upper().isin(list("ACGT")).all() or not variants["ALT"].str.upper().isin(list("ACGT")).all():
		raise ValueError("Only biallelic A/C/G/T SNPs are supported by this preparation workflow.")
	if (variants["REF"].str.upper() == variants["ALT"].str.upper()).any():
		raise ValueError("REF and ALT must differ.")
	labels = data[traits].apply(pd.to_numeric, errors="raise").to_numpy(dtype=np.float64)
	if np.isinf(labels).any():
		raise ValueError("Phenotypes must be finite or missing, never +/- infinity.")
	for column, trait in enumerate(traits):
		observed = labels[np.isfinite(labels[:, column]), column]
		if SUPPORTED_TRAITS[trait] == "binary" and not np.isin(observed, [0, 1]).all():
			raise ValueError("t2dm must be explicitly defined baseline case/control status coded 0/1/NA. Incident events are not substituted.")
		if require_training:
			for split in ("train", "validation"):
				values = labels[(data["split"].to_numpy() == split) & np.isfinite(labels[:, column]), column]
				if len(values) < 2:
					raise ValueError(f"{trait} needs at least two observed phenotypes in {split}.")
				if SUPPORTED_TRAITS[trait] == "binary" and len(np.unique(values)) < 2:
					raise ValueError(f"t2dm requires both cases and controls in {split}.")
	# This is an inexpensive format check; each actual batch is checked again.
	sample_rows = np.unique(np.linspace(0, len(data) - 1, min(24, len(data)), dtype=int))
	for row in sample_rows:
		validate_genotypes(genotypes[row])
	variant_identity = variants[variant_columns].to_csv(sep="\t", index=False)
	variant_digest = hashlib.sha256(variant_identity.encode()).hexdigest()
	data_digest = hashlib.sha256(data.to_csv(sep="\t", index=False).encode()).hexdigest()
	genotype_stat = Path(args.genotypes).stat()
	identity = {
		"data_sha256": data_digest, "variants_sha256": variant_digest,
		"data_file_sha256": file_sha256(args.data), "variants_file_sha256": file_sha256(args.variants),
		"genotypes_size": int(genotype_stat.st_size), "genotypes_mtime_ns": int(genotype_stat.st_mtime_ns),
		"shape": list(genotypes.shape), "genotypes_dtype": str(genotypes.dtype),
	}
	final_stamps = [(path.stat().st_size, path.stat().st_mtime_ns, path.stat().st_ino) for path in input_paths]
	if initial_stamps != final_stamps:
		raise RuntimeError("An input file changed while it was being read; retry with a stable prepared cache.")
	return data, labels, identity, variants[variant_columns]


def covariate_design(data, columns, training_rows=None, transform=None):
	missing = [column for column in columns if column not in data]
	if missing:
		raise ValueError(f"Covariate columns not found: {missing}")
	values = data[columns].apply(pd.to_numeric, errors="raise").to_numpy(dtype=np.float64)
	if np.isinf(values).any():
		raise ValueError("Covariates must be finite or missing, not infinity.")
	if transform is None:
		means, scales = [], []
		for column, name in enumerate(columns):
			observed = values[training_rows, column]
			observed = observed[np.isfinite(observed)]
			if not len(observed):
				raise ValueError(f"No observed training values for covariate {name}.")
			means.append(float(observed.mean()))
			scale = float(observed.std())
			scales.append(scale if scale > 1e-12 else 1.0)
		transform = {"columns": columns, "means": means, "scales": scales}
	means, scales = np.asarray(transform["means"]), np.asarray(transform["scales"])
	if columns:
		values = np.where(np.isfinite(values), values, means)
		values = (values - means) / scales
	design = np.column_stack([np.ones(len(data)), values])
	return design, transform


def fit_binary_baseline(design, phenotype, ridge=1e-4):
	# A small fixed penalty stabilizes covariate-only fits; the intercept is unpenalized.
	initial = np.zeros(design.shape[1])
	prevalence = float(phenotype.mean())
	initial[0] = math.log(prevalence / (1 - prevalence))

	def objective(coefficients):
		logits = design @ coefficients
		penalty = coefficients.copy()
		penalty[0] = 0
		loss = np.mean(np.logaddexp(0, logits) - phenotype * logits) + ridge * np.sum(penalty ** 2) / 2
		gradient = design.T @ (expit(logits) - phenotype) / len(phenotype) + ridge * penalty
		return float(loss), gradient

	result = minimize(objective, initial, jac=True, method="L-BFGS-B", options={"maxiter": 1000, "ftol": 1e-12, "gtol": 1e-7})
	if not result.success or not np.isfinite(result.x).all():
		raise ValueError(f"Covariate-only logistic regression did not converge: {result.message}")
	return result.x


def fit_baselines(data, labels, traits, columns):
	training_rows = np.flatnonzero(data["split"].to_numpy() == "train")
	design, transform = covariate_design(data, columns, training_rows=training_rows)
	baselines, baseline_logits = np.zeros_like(labels), np.zeros_like(labels)
	prepared_labels = np.full_like(labels, np.nan)
	task_parameters = {}
	for column, trait in enumerate(traits):
		fit_rows = training_rows[np.isfinite(labels[training_rows, column])]
		phenotype = labels[fit_rows, column]
		if SUPPORTED_TRAITS[trait] == "continuous":
			coefficients = np.linalg.lstsq(design[fit_rows], phenotype, rcond=None)[0]
			baseline = design @ coefficients
			scale = float(np.std(phenotype - baseline[fit_rows], ddof=0))
			if not np.isfinite(scale) or scale < 1e-12:
				raise ValueError(f"{trait} has no residual training variance after covariate adjustment.")
			baselines[:, column] = baseline
			prepared_labels[:, column] = (labels[:, column] - baseline) / scale
		else:
			coefficients = fit_binary_baseline(design[fit_rows], phenotype)
			baseline_logits[:, column] = design @ coefficients
			baselines[:, column] = expit(baseline_logits[:, column])
			prepared_labels[:, column] = labels[:, column]
			scale = 1.0
		task_parameters[trait] = {
			"type": SUPPORTED_TRAITS[trait], "coefficients": coefficients.tolist(),
			"residual_scale": scale, "n_training": int(len(fit_rows)),
			"logistic_ridge": 1e-4 if SUPPORTED_TRAITS[trait] == "binary" else None,
		}
	return {"covariates": transform, "tasks": task_parameters}, prepared_labels, baselines, baseline_logits


def apply_baselines(data, traits, preprocessing):
	transform = preprocessing["covariates"]
	design, _ = covariate_design(data, transform["columns"], transform=transform)
	baselines = np.zeros((len(data), len(traits)), dtype=np.float64)
	baseline_logits = np.zeros_like(baselines)
	for column, trait in enumerate(traits):
		linear_prediction = design @ np.asarray(preprocessing["tasks"][trait]["coefficients"])
		if SUPPORTED_TRAITS[trait] == "continuous":
			baselines[:, column] = linear_prediction
		else:
			baseline_logits[:, column] = linear_prediction
			baselines[:, column] = expit(linear_prediction)
	return baselines, baseline_logits


# 🚩 Official source loading; limited import compatibility without editing upstream
def load_official_model(upstream_dir, attention):
	upstream_dir = Path(upstream_dir).resolve()
	source_directory = upstream_dir / "src"
	paths = {name: source_directory / f"{name}.py" for name in ("utils", "modules", "model")}
	if not all(path.is_file() for path in paths.values()):
		raise ValueError("--upstream-dir must contain the official src/model.py, src/modules.py, and src/utils.py.")
	compatibility = []
	if attention == "neighborhood":
		try:
			natten = importlib.import_module("natten")
		except ImportError as exc:
			raise RuntimeError("Official neighborhood attention needs NATTEN. Use the documented torch 2.6.0 + natten 0.17.5 environment; no global-attention fallback is used.") from exc
		if str(getattr(natten, "__version__", "")).split("+")[0] != "0.17.5":
			raise RuntimeError("This adapter pins NATTEN 0.17.5 for the upstream use_fused_na/is_fna_enabled API; install the matching PyTorch/CUDA wheel.")
		if not all(hasattr(natten, name) for name in ("NeighborhoodAttention1D", "use_fused_na", "is_fna_enabled")):
			raise RuntimeError("The installed NATTEN build does not expose the API used by official PRSformer.")
	previous_modules = {name: sys.modules.get(name) for name in paths}
	loaded = {}
	try:
		for name, path in paths.items():
			tree = ast.parse(path.read_text(), filename=str(path))
			if name == "utils":
				lightning_imports = [node for node in tree.body if isinstance(node, ast.Import) and any(alias.name == "pytorch_lightning" for alias in node.names)]
				if lightning_imports:
					aliases = {alias.asname or alias.name for node in lightning_imports for alias in node.names if alias.name == "pytorch_lightning"}
					if any(isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load) and node.id in aliases for node in ast.walk(tree)):
						raise RuntimeError("Upstream now uses pytorch_lightning actively; review this import adapter before using that revision.")
					for node in lightning_imports:
						node.names = [alias for alias in node.names if alias.name != "pytorch_lightning"]
					tree.body = [node for node in tree.body if not isinstance(node, ast.Import) or node.names]
					compatibility.append("Skipped the unused pytorch_lightning import in utils.py; all functions are unchanged.")
			if name == "modules" and attention == "global":
				# These names occur only inside the unused neighborhood-attention class.
				tree.body = [node for node in tree.body if not (
					(isinstance(node, ast.ImportFrom) and node.module == "natten") or
					(isinstance(node, ast.Import) and any(alias.name == "natten" for alias in node.names))
				)]
				compatibility.append("Explicit small global-attention mode skips NATTEN imports and uses the upstream PyTorch attention branch.")
			module = types.ModuleType(f"grid_prsformer_official_{name}")
			module.__file__ = str(path)
			module.__package__ = ""
			sys.modules[name] = module
			sys.modules[module.__name__] = module
			exec(compile(ast.fix_missing_locations(tree), str(path), "exec"), module.__dict__)
			loaded[name] = module
	finally:
		for name, original in previous_modules.items():
			if original is None:
				sys.modules.pop(name, None)
			else:
				sys.modules[name] = original
	model_class = getattr(loaded["model"], "g2p_transformer_ExplicitNaNDose2", None)
	if model_class is None:
		raise RuntimeError("The official g2p_transformer_ExplicitNaNDose2 model class was not found.")
	commit = subprocess.run(["git", "-C", str(upstream_dir), "rev-parse", "HEAD"], capture_output=True, text=True, check=False)
	metadata = {
		"repository": "https://github.com/23andMe/PRSformer",
		"class": "g2p_transformer_ExplicitNaNDose2",
		"commit": commit.stdout.strip() if commit.returncode == 0 else None,
		"source_sha256": {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in paths.items()},
		"import_compatibility": compatibility,
	}
	return model_class, metadata


def architecture_from_args(args, length, tasks):
	dilation = [int(value) for value in csv_list(args.dilation)]
	if len(dilation) == 1:
		dilation *= args.layers
	if len(dilation) != args.layers or any(value < 1 for value in dilation):
		raise ValueError("--dilation needs one positive integer, or one positive integer per layer.")
	if min(length, tasks, args.embed_dim, args.heads, args.layers, args.ff_dim) < 1 or args.embed_dim % args.heads:
		raise ValueError("Model dimensions must be positive; --embed-dim must be divisible by --heads.")
	if args.attention == "global" and length > 4096:
		raise ValueError("Global attention is only an explicit small-data check (at most 4096 SNPs), not a genome-scale fallback.")
	if args.attention == "neighborhood" and (args.kernel_size < 2 or args.kernel_size * max(dilation) > length):
		raise ValueError("Neighborhood attention requires kernel_size >= 2 and kernel_size * max(dilation) <= the SNP count.")
	return {
		"seq_len": int(length), "embed_dim": args.embed_dim, "num_heads": args.heads,
		"dim_feedforward": args.ff_dim, "num_layers": args.layers, "num_covars": 0,
		"num_phenos": int(tasks), "kernel_size": args.kernel_size if args.attention == "neighborhood" else None,
		"dilation": dilation, "dlm_reprs": None, "weight_the_loss": False,
		"use_snp_annots": False, "snp_indices": None,
	}


def parameter_estimate(architecture):
	length, dimension = architecture["seq_len"], architecture["embed_dim"]
	feedforward, tasks = architecture["dim_feedforward"], architecture["num_phenos"]
	per_layer = 4 * dimension ** 2 + 2 * dimension * feedforward + 9 * dimension + feedforward
	return int(length * dimension * (2 + tasks) + architecture["num_layers"] * per_layer + 2 * dimension + tasks)


def import_torch():
	try:
		return importlib.import_module("torch")
	except ImportError as exc:
		raise RuntimeError("PyTorch is missing. Use the separate PRSformer environment documented with 3.prsformer.sh.") from exc


def runtime(args, attention):
	torch = import_torch()
	torch.set_num_threads(args.threads)
	device = torch.device(args.device)
	if device.type not in ("cuda", "cpu"):
		raise ValueError("This adapter supports cuda or cpu devices.")
	if device.type == "cuda" and not torch.cuda.is_available():
		raise RuntimeError("CUDA is unavailable; genome-scale PRSformer training needs a compatible NVIDIA GPU.")
	if device.type == "cuda":
		torch.cuda.set_device(device)
	if attention == "neighborhood" and device.type != "cuda":
		raise RuntimeError("The upstream fused neighborhood attention path requires CUDA. CPU is supported only for an explicit small global-attention check.")
	if args.amp == "bf16" and device.type == "cuda" and not torch.cuda.is_bf16_supported():
		raise RuntimeError("This GPU does not support BF16; select --amp fp16 or off.")
	if device.type == "cpu" and args.amp != "off":
		raise ValueError("For CPU small-data checks set --amp off explicitly.")
	dtype = {"fp16": torch.float16, "bf16": torch.bfloat16, "off": torch.float32}[args.amp]
	return torch, device, dtype


def checkpoint_layers(model, torch):
	from torch.utils.checkpoint import checkpoint
	for layer in model.transformer_encoder:
		original_forward = layer.forward

		def forward(value=None, *, x=None, _layer=layer, _forward=original_forward):
			inputs = x if x is not None else value
			if _layer.training and torch.is_grad_enabled():
				return checkpoint(_forward, inputs, use_reentrant=False)
			return _forward(inputs)

		layer.forward = forward


def model_smoke(model_class, architecture, torch, device, dtype, args):
	probe = dict(architecture)
	probe["seq_len"] = max(8, (architecture["kernel_size"] or 8) * max(architecture["dilation"]))
	model = model_class(**probe).to(device)
	if args.gradient_checkpointing:
		checkpoint_layers(model, torch)
	genotypes = torch.randint(0, 3, (2, probe["seq_len"]), device=device).float()
	genotypes[0, 0] = -1
	with torch.autocast(device_type=device.type, enabled=args.amp != "off", dtype=dtype):
		output = model(genotypes)
		loss = output.float().square().mean()
	if output.shape != (2, architecture["num_phenos"]) or not torch.isfinite(loss):
		raise RuntimeError("The official model failed its shape/finite-value preflight.")
	loss.backward()
	if not all(parameter.grad is None or torch.isfinite(parameter.grad).all() for parameter in model.parameters()):
		raise RuntimeError("The official model produced nonfinite gradients in the device preflight.")
	del model, genotypes, output, loss
	if device.type == "cuda":
		torch.cuda.empty_cache()


# 🚩 Bounded-memory batches, masked mixed-task losses and independent evaluation
class GenotypeRows:
	def __init__(self, filename, rows):
		self.filename = str(filename)
		self.rows = np.asarray(rows, dtype=np.int64)
		self.matrix = None

	def __len__(self):
		return len(self.rows)

	def __getitem__(self, index):
		if self.matrix is None:
			self.matrix = np.load(self.filename, mmap_mode="r", allow_pickle=False)
		row = int(self.rows[index])
		values = np.array(self.matrix[row], dtype=np.float32, copy=True)
		validate_genotypes(values)
		return values, row

	def __getstate__(self):
		state = dict(self.__dict__)
		state["matrix"] = None
		return state


def make_loader(args, rows, torch, shuffle=False, seed=None):
	generator = torch.Generator()
	generator.manual_seed(args.seed if seed is None else seed)
	return torch.utils.data.DataLoader(
		GenotypeRows(args.genotypes, rows), batch_size=args.batch_size, shuffle=shuffle,
		num_workers=args.workers, pin_memory=args.device.startswith("cuda"),
		generator=generator, drop_last=False,
	)


def loss_matrix(output, row_indices, prepared_labels, baseline_logits, traits, torch, device):
	target_array = prepared_labels[row_indices]
	mask = torch.as_tensor(np.isfinite(target_array), device=device)
	target = torch.as_tensor(np.nan_to_num(target_array, nan=0.0), device=device, dtype=torch.float32)
	offsets = torch.as_tensor(baseline_logits[row_indices], device=device, dtype=torch.float32)
	columns = []
	for column, trait in enumerate(traits):
		if SUPPORTED_TRAITS[trait] == "continuous":
			loss = (output[:, column].float() - target[:, column]) ** 2
		else:
			loss = torch.nn.functional.binary_cross_entropy_with_logits(
				output[:, column].float() + offsets[:, column], target[:, column], reduction="none",
			)
		columns.append(loss)
	return torch.stack(columns, dim=1) * mask, mask


def evaluate_loss(model, rows, args, prepared_labels, baseline_logits, traits, torch, device, dtype):
	model.eval()
	totals, counts = np.zeros(len(traits)), np.zeros(len(traits), dtype=np.int64)
	with torch.no_grad():
		for genotypes, row_indices in make_loader(args, rows, torch):
			genotypes = genotypes.to(device, non_blocking=True)
			row_indices = row_indices.numpy()
			with torch.autocast(device_type=device.type, enabled=args.amp != "off", dtype=dtype):
				output = model(genotypes)
				losses, mask = loss_matrix(output, row_indices, prepared_labels, baseline_logits, traits, torch, device)
			if not torch.isfinite(losses).all():
				raise RuntimeError("Nonfinite validation loss; no model is selected using invalid values.")
			totals += losses.sum(dim=0).cpu().numpy()
			counts += mask.sum(dim=0).cpu().numpy()
	return np.divide(totals, counts, out=np.full_like(totals, np.nan), where=counts > 0), counts


def predict_outputs(model, rows, args, torch, device, dtype, tasks):
	model.eval()
	outputs = np.empty((len(rows), tasks), dtype=np.float64)
	position = 0
	with torch.no_grad():
		for genotypes, _ in make_loader(args, rows, torch):
			with torch.autocast(device_type=device.type, enabled=args.amp != "off", dtype=dtype):
				prediction = model(genotypes.to(device, non_blocking=True))
			if not torch.isfinite(prediction).all():
				raise RuntimeError("The model produced nonfinite test predictions.")
			count = len(prediction)
			outputs[position:position + count] = prediction.float().cpu().numpy()
			position += count
	return outputs


def correlation_squared(first, second):
	if len(first) < 2 or np.std(first) < 1e-12 or np.std(second) < 1e-12:
		return np.nan
	return float(np.corrcoef(first, second)[0, 1] ** 2)


def auc(phenotype, prediction):
	positive = phenotype == 1
	n_positive, n_negative = int(positive.sum()), int((~positive).sum())
	if not n_positive or not n_negative:
		return np.nan
	ranks = rankdata(prediction, method="average")
	return float((ranks[positive].sum() - n_positive * (n_positive + 1) / 2) / (n_positive * n_negative))


def average_precision(phenotype, prediction):
	if not np.any(phenotype == 1) or not np.any(phenotype == 0):
		return np.nan
	order = np.argsort(-prediction, kind="mergesort")
	scores, labels = prediction[order], phenotype[order]
	endpoints = np.r_[np.flatnonzero(np.diff(scores) != 0), len(scores) - 1]
	true_positives = np.cumsum(labels)[endpoints]
	precision = true_positives / (endpoints + 1)
	recall = true_positives / labels.sum()
	return float(np.sum(np.diff(np.r_[0, recall]) * precision))


def metrics_for(phenotype, prediction, baseline, genetic_score, trait):
	if not len(phenotype):
		return {}
	if SUPPORTED_TRAITS[trait] == "continuous":
		baseline_error = np.sum((phenotype - baseline) ** 2)
		full_error = np.sum((phenotype - prediction) ** 2)
		total_variance = np.sum((phenotype - phenotype.mean()) ** 2)
		return {
			"prediction_R2": correlation_squared(phenotype - baseline, genetic_score),
			"SSE_partial_R2": float(1 - full_error / baseline_error) if baseline_error > 0 else np.nan,
			"full_R2": float(1 - full_error / total_variance) if total_variance > 0 else np.nan,
			"baseline_R2": float(1 - baseline_error / total_variance) if total_variance > 0 else np.nan,
			"full_RMSE": float(np.sqrt(full_error / len(phenotype))),
			"baseline_RMSE": float(np.sqrt(baseline_error / len(phenotype))),
			"prediction_bias": float(np.mean(prediction - phenotype)),
		}
	full_auc, baseline_auc = auc(phenotype, prediction), auc(phenotype, baseline)
	clipped = np.clip(prediction, 1e-12, 1 - 1e-12)
	return {
		"AUC": full_auc, "baseline_AUC": baseline_auc, "delta_AUC": full_auc - baseline_auc,
		"AUPRC": average_precision(phenotype, prediction), "baseline_AUPRC": average_precision(phenotype, baseline),
		"Brier": float(np.mean((phenotype - prediction) ** 2)),
		"baseline_Brier": float(np.mean((phenotype - baseline) ** 2)),
		"log_loss": float(-np.mean(phenotype * np.log(clipped) + (1 - phenotype) * np.log1p(-clipped))),
		"prevalence": float(phenotype.mean()), "cases": int(phenotype.sum()),
	}


def result_tables(data, labels, rows, raw_outputs, baselines, baseline_logits, traits, preprocessing):
	prediction_tables, metric_records = [], []
	for column, trait in enumerate(traits):
		baseline = baselines[rows, column]
		if SUPPORTED_TRAITS[trait] == "continuous":
			genetic = raw_outputs[:, column] * preprocessing["tasks"][trait]["residual_scale"]
			prediction = baseline + genetic
		else:
			genetic = raw_outputs[:, column]
			prediction = expit(baseline_logits[rows, column] + genetic)
		table = data.iloc[rows][["eid", "target", "split"]].copy()
		table["trait"] = trait
		table["y"] = labels[rows, column]
		table["baseline"], table["prediction"], table["genetic_score"] = baseline, prediction, genetic
		prediction_tables.append(table)
		groups = [("ALL", np.ones(len(table), dtype=bool))]
		groups += [(str(target), table["target"].to_numpy() == target) for target in sorted(table["target"].unique()) if target != "ALL"]
		for target, membership in groups:
			observed = membership & np.isfinite(table["y"].to_numpy())
			values = metrics_for(table["y"].to_numpy()[observed], prediction[observed], baseline[observed], genetic[observed], trait)
			for metric, value in values.items():
				metric_records.append({"split": "test", "target": target, "trait": trait, "n": int(observed.sum()), "metric": metric, "value": value})
	return pd.concat(prediction_tables, ignore_index=True), pd.DataFrame(metric_records)


# 🚩 Atomic checkpoints and test-only publication staging
def atomic_torch_save(torch, payload, filename):
	filename = Path(filename)
	staging = filename.with_name(filename.name + f".tmp.{os.getpid()}")
	try:
		with staging.open("wb") as handle:
			torch.save(payload, handle)
			handle.flush()
			os.fsync(handle.fileno())
		os.replace(staging, filename)
	finally:
		staging.unlink(missing_ok=True)


def atomic_table(table, filename):
	filename = Path(filename)
	staging = filename.with_name(filename.name + f".tmp.{os.getpid()}")
	try:
		table.to_csv(staging, sep="\t", index=False, compression="gzip" if str(filename).endswith(".gz") else None)
		os.replace(staging, filename)
	finally:
		staging.unlink(missing_ok=True)


def stage_evaluation(args, identity, preparation, checkpoint):
	# This completion record is temporary. Publication checks it and exports no JSON.
	genotype_stat = Path(args.genotypes).stat()
	if genotype_stat.st_size != identity["genotypes_size"] or genotype_stat.st_mtime_ns != identity["genotypes_mtime_ns"]:
		raise RuntimeError("The genotype cache changed while the model was running; outputs are not marked complete.")
	if file_sha256(args.data) != identity["data_file_sha256"] or file_sha256(args.variants) != identity["variants_file_sha256"]:
		raise RuntimeError("The input cohort or SNP manifest changed while the model was running; outputs are not marked complete.")
	output_directory = Path(args.out_dir)
	model_stat = (output_directory / "model.pt").stat()
	metadata = {
		"format_version": 1,
		"preparation_signature": preparation.get("signature") if preparation else None,
		"inputs": {
			"data_sha256": identity["data_file_sha256"], "variants_sha256": identity["variants_file_sha256"],
			"canonical_data_sha256": identity["data_sha256"], "canonical_variants_sha256": identity["variants_sha256"],
			"genotypes": {"size": identity["genotypes_size"], "mtime_ns": identity["genotypes_mtime_ns"], "shape": identity["shape"], "dtype": identity["genotypes_dtype"]},
		},
		"outputs_sha256": {name: file_sha256(output_directory / name) for name in ("metrics.tsv", "test_predictions.tsv.gz", "history.tsv")},
		"model": {
			"file": "model.pt", "size": int(model_stat.st_size), "mtime_ns": int(model_stat.st_mtime_ns),
			"training_signature": checkpoint["signature"], "upstream_source_sha256": checkpoint["upstream"]["source_sha256"],
		},
	}
	filename = output_directory / "evaluation.json"
	staging = filename.with_name(filename.name + f".tmp.{os.getpid()}")
	try:
		staging.write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n")
		os.replace(staging, filename)
	finally:
		staging.unlink(missing_ok=True)


def cpu_state(model):
	return {name: value.detach().cpu() for name, value in model.state_dict().items()}


@contextlib.contextmanager
def output_lock(directory):
	import fcntl
	key = hashlib.sha256(str(Path(directory).resolve()).encode()).hexdigest()[:20]
	with open(f"/tmp/grid-prsformer-{key}.lock", "a+") as handle:
		try:
			fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
		except BlockingIOError as exc:
			raise RuntimeError("Another PRSformer process is writing this output directory.") from exc
		yield


def save_history(history, filename, best_epoch):
	table = pd.DataFrame(history)
	table["selected"] = table["epoch"] == best_epoch
	atomic_table(table, filename)


def run_training(args):
	traits = csv_list(args.traits) or list(SUPPORTED_TRAITS)
	if len(set(traits)) != len(traits) or not set(traits).issubset(SUPPORTED_TRAITS):
		raise ValueError("--traits must be a nonrepeating subset of height,ldl,t2dm.")
	covariates = [] if args.covariates.strip().lower() == "none" else csv_list(args.covariates)
	if len(set(covariates)) != len(covariates):
		raise ValueError("Covariate names must not repeat.")
	if set(covariates) & set(traits):
		raise ValueError("A modeled phenotype cannot also be supplied as a covariate.")
	data, labels, identity, variants = read_inputs(args, traits)
	preparation_path = Path(args.genotypes).resolve().parent / "prepare.json"
	preparation = json.loads(preparation_path.read_text()) if preparation_path.is_file() else None
	architecture = architecture_from_args(args, identity["shape"][1], len(traits))
	preprocessing, prepared_labels, baselines, baseline_logits = fit_baselines(data, labels, traits, covariates)
	torch, device, dtype = runtime(args, args.attention)
	model_class, upstream = load_official_model(args.upstream_dir, args.attention)
	model_smoke(model_class, architecture, torch, device, dtype, args)
	parameter_count = parameter_estimate(architecture)
	summary = {
		"samples": len(data), "variants": identity["shape"][1], "traits": traits,
		"split_counts": {str(key): int(value) for key, value in data["split"].value_counts().items()},
		"observed_training": {trait: preprocessing["tasks"][trait]["n_training"] for trait in traits},
		"parameters": parameter_count, "fp32_parameters_GiB": parameter_count * 4 / 2 ** 30,
		"adam_parameters_gradients_states_minimum_GiB": parameter_count * 16 / 2 ** 30,
		"one_float32_activation_GiB": args.batch_size * identity["shape"][1] * args.embed_dim * 4 / 2 ** 30,
		"device": str(device), "attention": args.attention, "official_forward_backward_preflight": "passed",
		"note": "Activation/workspace memory is additional; parameter estimates do not guarantee a model will fit.",
	}
	print(json.dumps(summary, indent=2), flush=True)
	if args.check:
		return
	if device.type == "cuda":
		free_memory, _ = torch.cuda.mem_get_info(device)
		if parameter_count * 16 > free_memory * 0.8:
			raise RuntimeError("Model weights/gradients/Adam states alone would consume over 80% of free GPU memory; reduce SNPs/model dimensions before running.")
	output_directory = Path(args.out_dir)
	output_directory.mkdir(parents=True, exist_ok=True)
	best_file, training_file = output_directory / "model.pt", output_directory / "training.pt"
	if args.resume:
		if not training_file.is_file() or not best_file.is_file():
			raise ValueError("--resume needs both training.pt and model.pt from an interrupted run.")
	elif best_file.exists() or training_file.exists():
		raise ValueError("Output already contains a model. Use --resume for an interrupted run or choose a new output directory.")
	(output_directory / "evaluation.json").unlink(missing_ok=True)
	train_rows = np.flatnonzero((data["split"].to_numpy() == "train") & np.isfinite(labels).any(axis=1))
	validation_rows = np.flatnonzero(data["split"].to_numpy() == "validation")
	test_rows = np.flatnonzero(data["split"].to_numpy() == "test")
	task_counts = np.isfinite(prepared_labels[train_rows]).sum(axis=0)
	task_adjustment = torch.as_tensor(len(train_rows) / task_counts, device=device, dtype=torch.float32)
	# Uniform subject sampling plus this adjustment estimates equal mean loss per task.
	training_options = {
		name: getattr(args, name) for name in (
			"epochs", "patience", "min_delta", "lr", "weight_decay", "batch_size", "accumulation_steps",
			"clip_grad", "amp", "seed", "gradient_checkpointing",
		)
	}
	adapter_logic = ast.dump(ast.parse(Path(__file__).read_text()), include_attributes=False)
	adapter_digest = hashlib.sha256(adapter_logic.encode()).hexdigest()
	signature = digest_json({
		"identity": identity, "architecture": architecture, "traits": traits, "covariates": covariates,
		"training": training_options, "sources": upstream["source_sha256"], "adapter_logic_sha256": adapter_digest,
		"preprocessing": preprocessing, "preparation": preparation,
	})
	random.seed(args.seed)
	np.random.seed(args.seed)
	torch.manual_seed(args.seed)
	if device.type == "cuda":
		torch.cuda.manual_seed_all(args.seed)
	model = model_class(**architecture).to(device)
	if args.gradient_checkpointing:
		checkpoint_layers(model, torch)
	if sum(parameter.numel() for parameter in model.parameters()) != parameter_count:
		raise RuntimeError("Upstream parameter count changed; review the architecture before proceeding.")
	optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, betas=(0.9, 0.999), weight_decay=args.weight_decay)
	scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs, eta_min=args.lr * 0.02)
	scaler = torch.amp.GradScaler(device.type, enabled=device.type == "cuda" and args.amp == "fp16")
	start_epoch, best_epoch, best_loss, stale_epochs, history = 1, 0, float("inf"), 0, []
	metadata = {
		"format_version": 1, "traits": traits, "architecture": architecture, "attention": args.attention,
		"upstream": upstream, "preprocessing": preprocessing, "training_identity": identity,
		"variants": variants.to_dict(orient="list"), "preparation": preparation,
		"adapter_logic_sha256": adapter_digest,
		"training_options": training_options, "signature": signature, "torch_version": str(torch.__version__),
		"adaptation": "Mixed continuous/binary tasks; training-only covariate baselines; standardized residual MSE and fixed-logistic-offset BCE; equal task mean losses; validation early stopping.",
		"variant_dosage": "Raw count of the ALT allele in the aligned variant manifest; -1 denotes missing.",
		"test_usage": "Held out from neural training, baseline fitting, preprocessing, early stopping, and model selection.",
	}
	if args.resume:
		state = torch.load(training_file, map_location="cpu", weights_only=True)
		if state.get("signature") != signature:
			raise ValueError("Resume input, split, source, preprocessing or model/training settings differ from the interrupted run.")
		model.load_state_dict(state["model_state"])
		optimizer.load_state_dict(state["optimizer_state"])
		scheduler.load_state_dict(state["scheduler_state"])
		scaler.load_state_dict(state["scaler_state"])
		start_epoch, best_epoch, best_loss = state["epoch"] + 1, state["best_epoch"], state["best_loss"]
		stale_epochs, history = state["stale_epochs"], state["history"]
		torch.set_rng_state(state["torch_rng"])
		if device.type == "cuda" and state["cuda_rng"]:
			torch.cuda.set_rng_state_all(state["cuda_rng"])
		del state
	for epoch in range(start_epoch, args.epochs + 1):
		if stale_epochs >= args.patience:
			break
		model.train()
		optimizer.zero_grad(set_to_none=True)
		started = time.monotonic()
		loss_totals, observed_counts = np.zeros(len(traits)), np.zeros(len(traits), dtype=np.int64)
		loader = make_loader(args, train_rows, torch, shuffle=True, seed=args.seed + epoch)
		learning_rate = float(optimizer.param_groups[0]["lr"])
		for batch, (genotypes, row_indices) in enumerate(loader):
			row_indices = row_indices.numpy()
			group_start = (batch // args.accumulation_steps) * args.accumulation_steps * args.batch_size
			effective_count = min(args.accumulation_steps * args.batch_size, len(train_rows) - group_start)
			with torch.autocast(device_type=device.type, enabled=args.amp != "off", dtype=dtype):
				output = model(genotypes.to(device, non_blocking=True))
				losses, mask = loss_matrix(output, row_indices, prepared_labels, baseline_logits, traits, torch, device)
				loss = (losses * task_adjustment).sum() / effective_count / len(traits)
			if not torch.isfinite(loss):
				raise RuntimeError("Nonfinite training loss. Reduce the learning rate or use FP32; an invalid model is not saved.")
			scaler.scale(loss).backward()
			loss_totals += losses.detach().sum(dim=0).cpu().numpy()
			observed_counts += mask.sum(dim=0).cpu().numpy()
			if (batch + 1) % args.accumulation_steps == 0 or batch + 1 == len(loader):
				scaler.unscale_(optimizer)
				# FP16 loss scaling can initially overflow; GradScaler records this,
				# skips the update and reduces its scale. FP32/BF16 failures are errors.
				torch.nn.utils.clip_grad_norm_(model.parameters(), args.clip_grad, error_if_nonfinite=not scaler.is_enabled())
				scaler.step(optimizer)
				scaler.update()
				optimizer.zero_grad(set_to_none=True)
		train_losses = loss_totals / observed_counts
		validation_losses, validation_counts = evaluate_loss(model, validation_rows, args, prepared_labels, baseline_logits, traits, torch, device, dtype)
		validation_loss = float(np.mean(validation_losses))
		seconds = float(time.monotonic() - started)
		for split, losses, counts in (("train", train_losses, observed_counts), ("validation", validation_losses, validation_counts)):
			for column, trait in enumerate(traits):
				history.append({"epoch": epoch, "split": split, "trait": trait, "loss": float(losses[column]), "n": int(counts[column]), "learning_rate": learning_rate, "seconds": seconds})
		if validation_loss < best_loss - args.min_delta:
			best_loss, best_epoch, stale_epochs = validation_loss, epoch, 0
			atomic_torch_save(torch, {**metadata, "model_state": cpu_state(model), "best_epoch": best_epoch, "best_validation_loss": best_loss}, best_file)
		else:
			stale_epochs += 1
		scheduler.step()
		save_history(history, output_directory / "history.tsv", best_epoch)
		atomic_torch_save(torch, {
			"signature": signature, "epoch": epoch, "best_epoch": best_epoch, "best_loss": best_loss,
			"stale_epochs": stale_epochs, "history": history, "model_state": cpu_state(model),
			"optimizer_state": optimizer.state_dict(), "scheduler_state": scheduler.state_dict(), "scaler_state": scaler.state_dict(),
			"torch_rng": torch.get_rng_state(), "cuda_rng": torch.cuda.get_rng_state_all() if device.type == "cuda" else [],
		}, training_file)
		print(f"Epoch {epoch}: train={float(np.mean(train_losses)):.6f}; validation={validation_loss:.6f}; best={best_epoch}; seconds={seconds:.1f}", flush=True)
	if not best_file.is_file():
		raise RuntimeError("Training did not produce a valid validation-selected model.")
	best = torch.load(best_file, map_location="cpu", weights_only=True)
	model.load_state_dict(best["model_state"])
	outputs = predict_outputs(model, test_rows, args, torch, device, dtype, len(traits))
	predictions, metrics = result_tables(data, labels, test_rows, outputs, baselines, baseline_logits, traits, preprocessing)
	atomic_table(predictions, output_directory / "test_predictions.tsv.gz")
	atomic_table(metrics, output_directory / "metrics.tsv")
	best["completed_epochs"] = max(record["epoch"] for record in history)
	best["test_evaluated"] = True
	best["history"] = history
	atomic_torch_save(torch, best, best_file)
	stage_evaluation(args, identity, preparation, best)
	training_file.unlink(missing_ok=True)
	print(f"Completed: selected epoch {best_epoch}; test-only predictions for {len(test_rows)} individuals.", flush=True)


def run_prediction(args):
	if not args.checkpoint:
		raise ValueError("predict requires --checkpoint model.pt.")
	checkpoint_path = Path(args.checkpoint).resolve()
	checkpoint_stat = checkpoint_path.stat()
	torch = import_torch()
	checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
	traits, architecture = checkpoint["traits"], checkpoint["architecture"]
	if args.traits and csv_list(args.traits) != traits:
		raise ValueError("Prediction traits and their order must match the checkpoint.")
	data, labels, identity, _ = read_inputs(args, traits, require_training=False)
	preparation_path = Path(args.genotypes).resolve().parent / "prepare.json"
	preparation = json.loads(preparation_path.read_text()) if preparation_path.is_file() else None
	if identity["variants_sha256"] != checkpoint["training_identity"]["variants_sha256"]:
		raise ValueError("The ordered SNP/REF/ALT manifest differs from model training; do not score misaligned genotypes.")
	torch, device, dtype = runtime(args, checkpoint["attention"])
	model_class, upstream = load_official_model(args.upstream_dir, checkpoint["attention"])
	if upstream["source_sha256"] != checkpoint["upstream"]["source_sha256"]:
		raise ValueError("Upstream source files changed since model training; restore the recorded revision before prediction.")
	model_smoke(model_class, architecture, torch, device, dtype, args)
	if args.check:
		print("Prediction input identities and official model forward/backward preflight passed.", flush=True)
		return
	model = model_class(**architecture).to(device)
	model.load_state_dict(checkpoint["model_state"])
	baselines, baseline_logits = apply_baselines(data, traits, checkpoint["preprocessing"])
	rows = np.flatnonzero(data["split"].to_numpy() == "test")
	outputs = predict_outputs(model, rows, args, torch, device, dtype, len(traits))
	predictions, metrics = result_tables(data, labels, rows, outputs, baselines, baseline_logits, traits, checkpoint["preprocessing"])
	output_directory = Path(args.out_dir).resolve()
	output_directory.mkdir(parents=True, exist_ok=True)
	(output_directory / "evaluation.json").unlink(missing_ok=True)
	atomic_table(predictions, output_directory / "test_predictions.tsv.gz")
	atomic_table(metrics, output_directory / "metrics.tsv")
	if not checkpoint.get("history"):
		raise ValueError("The checkpoint lacks training history; use a complete checkpoint produced by this adapter.")
	save_history(checkpoint["history"], output_directory / "history.tsv", checkpoint["best_epoch"])
	model_destination = output_directory / "model.pt"
	if checkpoint_path != model_destination:
		staging = model_destination.with_name(model_destination.name + f".tmp.{os.getpid()}")
		try:
			shutil.copy2(checkpoint_path, staging)
			os.replace(staging, model_destination)
		finally:
			staging.unlink(missing_ok=True)
	current_checkpoint_stat = checkpoint_path.stat()
	if current_checkpoint_stat.st_size != checkpoint_stat.st_size or current_checkpoint_stat.st_mtime_ns != checkpoint_stat.st_mtime_ns:
		raise RuntimeError("The source checkpoint changed during prediction; outputs are not marked complete.")
	stage_evaluation(args, identity, preparation, checkpoint)
	print(f"Scored {len(rows)} held-out individuals using the saved model and training-only preprocessing.", flush=True)


# 🚩 Command-line entry point
def main():
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("command", choices=["train", "predict"])
	parser.add_argument("--genotypes", required=True, help="Aligned sample-major .npy dosage memmap")
	parser.add_argument("--data", required=True, help="TSV(.gz): eid,target,split,phenotypes,numeric covariates")
	parser.add_argument("--variants", required=True, help="Ordered TSV(.gz): CHR,BP,SNP,REF,ALT")
	parser.add_argument("--upstream-dir", required=True, help="Official 23andMe/PRSformer checkout")
	parser.add_argument("--out-dir", required=True, help="Private staging directory; publish temporary TSV files using 3.prsformer.sh")
	parser.add_argument("--traits", help="Default height,ldl,t2dm; prediction uses the checkpoint task order")
	parser.add_argument("--covariates", default="", help="Comma-separated numeric covariates; transformations fit only on training samples")
	parser.add_argument("--device", default="cuda")
	parser.add_argument("--attention", choices=["neighborhood", "global"], default="neighborhood")
	parser.add_argument("--embed-dim", type=int, default=64)
	parser.add_argument("--heads", type=int, default=4)
	parser.add_argument("--layers", type=int, default=2)
	parser.add_argument("--ff-dim", type=int, default=128)
	parser.add_argument("--kernel-size", type=int, default=385)
	parser.add_argument("--dilation", default="1")
	parser.add_argument("--batch-size", type=int, default=1)
	parser.add_argument("--accumulation-steps", type=int, default=64)
	parser.add_argument("--epochs", type=int, default=30)
	parser.add_argument("--patience", type=int, default=5)
	parser.add_argument("--min-delta", type=float, default=1e-4)
	parser.add_argument("--lr", type=float, default=5e-4)
	parser.add_argument("--weight-decay", type=float, default=0.05)
	parser.add_argument("--clip-grad", type=float, default=1.0)
	parser.add_argument("--amp", choices=["fp16", "bf16", "off"], default="fp16")
	parser.add_argument("--gradient-checkpointing", action=argparse.BooleanOptionalAction, default=True)
	parser.add_argument("--workers", type=int, default=0)
	parser.add_argument("--threads", type=int, default=4)
	parser.add_argument("--seed", type=int, default=12345)
	parser.add_argument("--check", action="store_true", help="Validate aligned inputs, baselines and a small actual official-model forward/backward pass")
	parser.add_argument("--resume", action="store_true", help="Resume an interrupted training run with identical settings")
	parser.add_argument("--checkpoint", help="Saved model.pt for predict")
	args = parser.parse_args()
	if min(args.batch_size, args.accumulation_steps, args.epochs, args.patience, args.threads) < 1 or args.workers < 0:
		parser.error("Batch/accumulation/epoch/patience/thread counts must be positive; workers must be nonnegative.")
	if not all(math.isfinite(value) for value in (args.lr, args.weight_decay, args.clip_grad, args.min_delta)) or args.lr <= 0 or args.weight_decay < 0 or args.clip_grad <= 0 or args.min_delta < 0:
		parser.error("Learning rate and gradient clip must be positive; weight decay and minimum improvement must be nonnegative.")
	with output_lock(args.out_dir):
		if args.command == "train":
			run_training(args)
		else:
			run_prediction(args)


if __name__ == "__main__":
	try:
		main()
	except (ValueError, RuntimeError, FileNotFoundError) as error:
		raise SystemExit(f"ERROR: {error}") from error
