import torch
import pickle
import numpy as np
from collections import Counter
import pandas as pd
import torch.nn.functional as F
from torch.autograd import Variable
import torch.nn as nn
from sklearn.metrics import recall_score, precision_score, f1_score
from torchvision import transforms


def ancestry_accuracy(prediction, target):
    b, l, c = prediction.shape
    prediction = prediction.reshape(b * l, c)
    target = target.reshape(b * l)

    prediction = prediction.max(dim=1)[1]
    accuracy = (prediction == target).sum()

    return accuracy / l



def filter_loci(batch):
    mixed_vcf = batch['mixed_vcf']
    mixed_labels = batch['mixed_labels']
    pos = batch['pos']
    ref_panel = batch['ref_panel']

    if isinstance(ref_panel, list) and all(isinstance(ref, dict) for ref in ref_panel):
        try:
            african_tensor = ref_panel[0][0].float()  
            den_tensor = ref_panel[0][1].float()      
            nean_tensor = ref_panel[0][2].float()     
        except Exception as e:
            print(f"Error when accessing ref_panel tensors: {e}")
            return batch  
    else:
        raise ValueError("Error")

    original_num_sites = mixed_vcf.shape[1]


    mask = torch.ones(mixed_vcf.shape[1], dtype=torch.bool, device=mixed_vcf.device)

    for idx in range(mixed_vcf.shape[1]):
        target_value = mixed_vcf[:, idx]
        african_value = african_tensor[:, idx].mean() if african_tensor.ndim > 1 else african_tensor[idx].float()
        den_value = den_tensor[:, idx].mean() if den_tensor.ndim > 1 else den_tensor[idx].float()
        nean_value = nean_tensor[:, idx].mean() if nean_tensor.ndim > 1 else nean_tensor[idx].float()
        if torch.all(target_value == african_value) and torch.all(target_value == den_value) and torch.all(target_value == nean_value):
            mask[idx] = False

    filtered_vcf = mixed_vcf[:, mask]
    filtered_labels = mixed_labels[:, mask]
    filtered_pos = pos[:, mask]
    filtered_num_sites = filtered_vcf.shape[1]

    batch['mixed_vcf'] = filtered_vcf
    batch['mixed_labels'] = filtered_labels
    batch['pos'] = filtered_pos

    for ref in ref_panel:
        try:
            ref[0] = ref[0][:, mask]  
            ref[1] = ref[1][:, mask]  
            ref[2] = ref[2][:, mask]  
        except Exception as e:
            print(f"Error when filtering ref_panel tensors: {e}")

    batch['ref_panel'] = ref_panel

    return batch


def ancestry_metrics(prediction, target, binary=False):
    prediction = prediction.permute(0,2,1)
    b, l, c = prediction.shape

    if binary:
        prediction = prediction.reshape(b*l)
        prediction = torch.floor(prediction + 0.5)
    else:
        prediction = prediction.reshape(b * l, c)
        prediction = prediction.max(dim=1)[1]
    target = target.reshape(b*l)

    
    accuracy = (prediction == target).sum()
    accuracy = accuracy / l /b
    target = target.cpu().numpy()
    prediction = prediction.cpu().numpy()
    recall = recall_score(target, prediction, average='macro', zero_division=0)
    precision = precision_score(target, prediction, average='macro',zero_division=0)
    f1 = f1_score(target, prediction, average='macro', zero_division=0)
    return accuracy, recall, precision, f1


def ancestry_metrics_label_based(prediction, target, binary=False):
    prediction = prediction.view(-1)
    target = target.view(-1)

    if prediction.ndimension() == 1 and prediction.shape[0] == target.shape[0]:
        accuracy = (prediction == target).sum().float() / prediction.size(0)

        target = target.cpu().numpy()
        prediction = prediction.cpu().numpy()

        recall = recall_score(target, prediction, average='macro', zero_division=0)
        precision = precision_score(target, prediction, average='macro', zero_division=0)
        f1 = f1_score(target, prediction, average='macro', zero_division=0)

        return accuracy, recall, precision, f1
    else:
        raise ValueError("Predictions and targets have incompatible shapes.")



def ancestry_metrics_bin(prediction, target, binary=False):
    prediction = prediction.permute(0,2,1)
    b, l, c = prediction.shape

    if binary:
        prediction = prediction.reshape(b*l)
        prediction = torch.floor(prediction + 0.5)
    else:
        prediction = prediction.reshape(b * l, c)
        prediction = prediction.max(dim=1)[1]
    target = target.reshape(b*l)

    
    accuracy = (prediction == target).sum()
    accuracy = accuracy / l /b
    target = target.cpu().numpy()
    prediction = prediction.cpu().numpy()
    recall = recall_score(target, prediction, zero_division=0)
    precision = precision_score(target, prediction,zero_division=0)
    f1 = f1_score(target, prediction, zero_division=0)
    return accuracy, recall, precision, f1

def ancestry_metrics_ad(output):
    b,l,s = output["predictions"].shape
    mse_neg = torch.mean(torch.pow(output["test"][output["train_indices"],0].reshape(-1,l*s) - output["test_predictions"][output['train_indices']].reshape(-1,l*s), 2), axis=1)
    mse_pos = torch.mean(torch.pow(output["test"][output["test_indices"],0].reshape(-1,l*s) - output["test_predictions"][output['test_indices']].reshape(-1,l*s), 2), axis=1)
    mean_neg = torch.mean(mse_neg)
    var_neg = torch.var(mse_neg)
    mean_pos = torch.mean(mse_pos)
    var_pos = torch.var(mse_pos)
    return mean_neg,var_neg,mean_pos,var_pos


class AverageMeter():
    def __init__(self):
        self.total = 0
        self.count = 0

    def reset(self):
        self.total = 0
        self.count = 0

    def update(self, value):
        self.total += value
        self.count += 1

    def get_average(self):
        return self.total / self.count


class ProgressSaver():

    def __init__(self, exp_dir):
        self.exp_dir = exp_dir
        self.progress = {
            "epoch": [],
            "train_loss": [],
            "val_loss": [],
            "val_acc": [],
            "val_recall": [],
            "val_precision": [],
            "val_f1": [],
            "time": [],
            "best_epoch": [],
            "best_val_f1": [],
            "best_val_loss": [],
            "lr": [],
            "iter": []
        }

    def update_epoch_progess(self, epoch_data):
        for key in epoch_data.keys():
            self.progress[key].append(epoch_data[key])

        with open("%s/progress.pckl" % self.exp_dir, "wb") as f:
            pickle.dump(self.progress, f)

    def load_progress(self):
        with open("%s/progress.pckl" % self.exp_dir, "rb") as f:
            self.progress = pickle.load(f)

    def get_resume_stats(self):
        return self.progress["best_epoch"][-1], self.progress["best_val_loss"][-1], self.progress["iter"][-1], self.progress["time"][-1]


class ReshapedCrossEntropyLoss(nn.Module):
    def __init__(self, loss):
        super(ReshapedCrossEntropyLoss, self).__init__()
        self.CELoss = nn.CrossEntropyLoss()
        self.BCELoss = nn.BCELoss()
        self.Focal = FocalLoss()
        self.mse = nn.L1Loss()
        self.loss = loss

    def forward(self, prediction, target, epoch):
        if self.loss == "MSE":
            loss = self.mse(prediction, target[:,0,:].unsqueeze(1)) #target[:,0,:].unsqueeze(1)
            return loss
        prediction = prediction.permute(0,2,1)
        bs, seq_len, n_classes = prediction.shape
        prediction = prediction.reshape(bs * seq_len, n_classes)

        target = target.reshape(bs * seq_len)
        if self.loss == "CE":
            loss = self.CELoss(prediction, target)
        elif self.loss == "BCE":
            target = target.reshape(bs * seq_len, 1).to(torch.float)
            loss = self.BCELoss(prediction, target)
        elif self.loss == "LDAM":
            cls_num_list = []
            cls_num_list.append((target==0).sum().item())
            cls_num_list.append((target==1).sum().item())
            for i in range(len(cls_num_list)):
                if cls_num_list[i] == 0:
                    cls_num_list[i] == 1
            loss = LDAMLoss(cls_num_list, weight=None).forward(prediction, target)
        elif self.loss == "VS":
            target = target.reshape(bs * seq_len, 1).to(torch.float)
            cls_num_list = []
            cls_num_list.append((target==0).sum().item())
            cls_num_list.append((target==1).sum().item())
            loss = VSLoss(cls_num_list).forward(prediction, target)
        else:
            loss = self.Focal(prediction, target)
        return loss


class LDAMLoss(nn.Module):
    
    def __init__(self, cls_num_list, max_m=0.5, weight=None, s=30):
        super(LDAMLoss, self).__init__()
        m_list = 1.0 / np.sqrt(np.sqrt(cls_num_list))
        m_list = m_list * (max_m / np.max(m_list))
        m_list = torch.cuda.FloatTensor(m_list)
        self.m_list = m_list
        assert s > 0
        self.s = s
        self.weight = weight

    def forward(self, x, target):
        target = target.type(torch.LongTensor).to(x.device)
        target = target.reshape(target.shape[0])
        index = torch.zeros_like(x, dtype=torch.uint8)
        index.scatter_(1, target.data.view(-1, 1), 1)
        
        index_float = index.type(torch.cuda.FloatTensor)
        batch_m = torch.matmul(self.m_list[None, :].to(index_float.device), index_float.transpose(0,1))
        batch_m = batch_m.view((-1, 1))
        x_m = x - batch_m
    
        output = torch.where(index, x_m, x)
        return F.cross_entropy(self.s*output, target, weight=self.weight)

class FocalLoss(nn.Module):
    def __init__(self, alpha=[0.25, 0.25, 0.25], gamma=2, logits=True, reduce=True):
        super(FocalLoss, self).__init__()
        self.alpha = torch.tensor(alpha)
        self.gamma = gamma
        self.logits = logits
        self.reduce = reduce

    def forward(self, inputs, targets):
        if self.logits:
            BCE_loss = F.cross_entropy(inputs, targets, reduction='none')
        else:
            BCE_loss = F.nll_loss(torch.log(inputs), targets, reduction='none')

        at =  self.alpha.to(targets.device)[targets]

        # gather the specific corresponding `-log(pt)` for each target class
        pt = torch.exp(-BCE_loss)
        F_loss = at * (1-pt)**self.gamma * BCE_loss

        if self.reduce:
            return torch.mean(F_loss)
        else:
            return F_loss


class VSLoss(nn.Module):

    def __init__(self, cls_num_list, gamma=0.2, tau=1.2, weight=None):
        super(VSLoss, self).__init__()

        cls_probs = [cls_num / sum(cls_num_list) for cls_num in cls_num_list]
        temp = (1.0 / np.array(cls_num_list)) ** gamma
        temp = temp / np.min(temp)

        iota_list = tau * np.log(cls_probs)
        Delta_list = temp

        self.iota_list = torch.cuda.FloatTensor(iota_list)
        self.Delta_list = torch.cuda.FloatTensor(Delta_list)
        self.weight = weight

    def forward(self, x, target):
        target = target.type(torch.LongTensor).to(x.device)
        target = target.reshape(target.shape[0])
        output = x / self.Delta_list + self.iota_list

        return F.cross_entropy(output, target, weight=self.weight)

def adjust_learning_rate(base_lr, lr_decay, optimizer, epoch):
    """Sets the learning rate to the initial LR decayed by 10 every lr_decay epochs"""
    lr = base_lr * (0.1 ** (epoch / lr_decay))
    for param_group in optimizer.param_groups:
        param_group['lr'] = lr

    return lr


class EncodeBinary:

    def __call__(self, inp):
        # 0 -> -1
        # 1 -> 1
        inp["mixed_vcf"] = inp["mixed_vcf"] * 2 - 1
        for anc in inp["ref_panel"]:
            inp["ref_panel"][anc] = inp["ref_panel"][anc] * 2 - 1

        return inp


def build_transforms(args):
    transforms_list = []

    transforms_list.append(EncodeBinary())

    transforms_list = transforms.Compose(transforms_list)

    return transforms_list


def to_device(item, device):
    item["mixed_vcf"] = item["mixed_vcf"].to(device)

    if "mixed_labels" in item.keys():
        item["mixed_labels"] = item["mixed_labels"].to(device)

    for i, panel in enumerate(item["ref_panel"]):
        for anc in panel.keys():
            item["ref_panel"][i][anc] = item["ref_panel"][i][anc].to(device)

    return item


def correct_max_indices(max_indices_batch, ref_panel_idx_batch):
    '''
    for each element of a batch, the dataloader samples randomly a set of founders in random order. For this reason,
    the argmax values output by the base model will represent different associations of founders, depending on how they have been
    sampled and ordered. By storing the sampling information during the data loading, we can then correct the argmax outputs
    into a shared meaning between batches and elements within the batch.
    '''

    for n in range(len(max_indices_batch)):

        max_indices = max_indices_batch[n]
        ref_panel_idx = ref_panel_idx_batch[n]
        max_indices_ordered = [None] * len(ref_panel_idx.keys())

        for i, c in enumerate(ref_panel_idx.keys()):
            max_indices_ordered[c] = max_indices[i]
        max_indices_ordered = torch.stack(max_indices_ordered)

        for i in range(max_indices.shape[0]):
            max_indices_ordered[i] = torch.take(torch.tensor(ref_panel_idx[i]), max_indices_ordered[i].cpu())

        max_indices_batch[n] = max_indices_ordered[:]

    return max_indices_batch


def compute_ibd(output):
    all_ibd = []
    for n in range(output['out_basemodel'].shape[0]):
        classes_basemodel = torch.argmax(output['out_basemodel'][n], dim=0)
        # classes_smoother = torch.argmax(output['out_smoother'][n], dim=0)
        ibd = torch.gather(output['max_indices'][n].t(), index=classes_basemodel.unsqueeze(1), dim=1)
        ibd = ibd.squeeze(1)

        all_ibd.append(ibd)

    all_ibd = torch.stack(all_ibd)

    return all_ibd


def get_meta_data(chm, model_pos, query_pos, n_wind, wind_size, gen_map_df=None):
    """
    from LAI-Net code
    Transforms the predictions on a window level to a .msp file format.
        - chm: chromosome number
        - model_pos: physical positions of the model input SNPs in basepair units
        - query_pos: physical positions of the query input SNPs in basepair units
        - n_wind: number of windows in model
        - wind_size: size of each window in the model
        - genetic_map_file: the input genetic map file
    """

    model_chm_len = len(model_pos)

    # chm
    chm_array = [chm] * n_wind

    # start and end pyshical positions
    if model_chm_len % wind_size == 0:
        spos_idx = np.arange(0, model_chm_len, wind_size)  # [:-1]
        epos_idx = np.concatenate([np.arange(0, model_chm_len, wind_size)[1:], np.array([model_chm_len])]) - 1
    else:
        spos_idx = np.arange(0, model_chm_len, wind_size)[:-1]
        epos_idx = np.concatenate([np.arange(0, model_chm_len, wind_size)[1:-1], np.array([model_chm_len])]) - 1

    spos = model_pos[spos_idx]
    epos = model_pos[epos_idx]

    sgpos = [1] * len(spos)
    egpos = [1] * len(epos)

    # number of query snps in interval
    wind_index = [min(n_wind - 1, np.where(q == sorted(np.concatenate([epos, [q]])))[0][0]) for q in query_pos]
    window_count = Counter(wind_index)
    n_snps = [window_count[w] for w in range(n_wind)]

    # print(len(chm_array), len(spos), len(epos), len(sgpos), len(egpos), len(n_snps))
    # Concat with prediction table
    meta_data = np.array([chm_array, spos, epos, sgpos, egpos, n_snps]).T
    meta_data_df = pd.DataFrame(meta_data)
    meta_data_df.columns = ["chm", "spos", "epos", "sgpos", "egpos", "n snps"]

    return meta_data_df


def write_msp_tsv(output_folder, meta_data, pred_labels, populations, query_samples, write_population_code=False):
    msp_data = np.concatenate([np.array(meta_data), pred_labels.T], axis=1).astype(str)

    with open(output_folder + "/predictions.msp.tsv", 'w') as f:
        if write_population_code:
            # first line (comment)
            f.write("#Subpopulation order/codes: ")
            f.write("\t".join([str(pop) + "=" + str(i) for i, pop in enumerate(populations)]) + "\n")
        # second line (comment/header)
        f.write("#" + "\t".join(meta_data.columns) + "\t")
        f.write("\t".join([str(s) for s in np.concatenate([[s + ".0", s + ".1"] for s in query_samples])]) + "\n")
        # rest of the lines (data)
        for l in range(msp_data.shape[0]):
            f.write("\t".join(msp_data[l, :]))
            f.write("\n")

    return


def msp_to_lai(msp_file, positions, lai_file=None):
    msp_df = pd.read_csv(msp_file, sep="\t", comment="#", header=None)
    data_window = np.array(msp_df.iloc[:, 6:])
    n_reps = msp_df.iloc[:, 5].to_numpy()
    assert np.sum(n_reps) == len(positions)
    data_snp = np.concatenate([np.repeat([row], repeats=n_reps[i], axis=0) for i, row in enumerate(data_window)])

    with open(msp_file) as f:
        first_line = f.readline()
        second_line = f.readline()

    header = second_line[:-1].split("\t")
    samples = header[6:]
    df = pd.DataFrame(data_snp, columns=samples, index=positions)

    if lai_file is not None:
        with open(lai_file, "w") as f:
            f.write(first_line)
        df.to_csv(lai_file, sep="\t", mode='a', index_label="position")




# Latest no overlap
def find_introgression_segments(
    df: pd.DataFrame,
    haplotype_columns: list,
    probabilities,
    Chr: int = 1,
    merge_distance: int = 0,
    max_snp_gap_threshold: int = 1_000_000,
    min_snps_per_segment: int = 2,
    mosaic_minority_threshold: float = 0.20
) -> pd.DataFrame:
    """
    Definitive, refactored function to identify and merge archaic introgression segments.
    This version uses a single-pass, ordered merging strategy to prevent any
    possibility of nested or overlapping segment outputs.
    """

    final_output_columns = [
        'chr', 'start_pos', 'end_pos', 'haplotype', 'label', 'snps',
        'prob', 'n_snps_label1', 'n_snps_label2'
    ]

    if not isinstance(df, pd.DataFrame) or df.empty or not all(c in df.columns for c in ['POS'] + haplotype_columns):
        return pd.DataFrame(columns=final_output_columns)

    pos_array = df['POS'].to_numpy(dtype=np.int64)
    num_snps_total = len(pos_array)

    try:
        prob_label1 = np.array(probabilities[0][1], dtype=np.float64)
        prob_label2 = np.array(probabilities[0][2], dtype=np.float64)
        if len(prob_label1) != num_snps_total or len(prob_label2) != num_snps_total:
            raise ValueError("Probability array lengths do not match SNP data length.")
    except Exception as e:
        raise ValueError(f"Error processing probabilities: {e}")

    all_haplotypes_final_segments = []

    for hap_col in haplotype_columns:
        hap_state_array = df[hap_col].to_numpy(dtype=np.int8)
        is_archaic_snp = np.isin(hap_state_array, [1, 2])

        if not np.any(is_archaic_snp):
            continue

        # --- 1. Find initial candidate blocks (non-overlapping) ---
        padded = np.concatenate(([False], is_archaic_snp, [False]))
        diffs = np.diff(padded.astype(np.int8))
        block_starts = np.where(diffs == 1)[0]
        block_ends = np.where(diffs == -1)[0] - 1

        candidate_blocks = []
        for s_idx, e_idx in zip(block_starts, block_ends):
            sub_block_start = s_idx
            positions_in_block = pos_array[s_idx:e_idx + 1]
            gaps = np.diff(positions_in_block)
            split_points = np.where(gaps > max_snp_gap_threshold)[0]
            
            for split_idx in split_points:
                sub_block_end = s_idx + split_idx
                if (sub_block_end - sub_block_start + 1) >= min_snps_per_segment:
                    candidate_blocks.append((sub_block_start, sub_block_end))
                sub_block_start = s_idx + split_idx + 1
            
            if (e_idx - sub_block_start + 1) >= min_snps_per_segment:
                candidate_blocks.append((sub_block_start, e_idx))
        
        if not candidate_blocks:
            continue

        # --- 2. Create a DataFrame of initial, classified, non-overlapping segments ---
        initial_segments = []
        for s_idx, e_idx in candidate_blocks:
            initial_segments.append({
                'start_pos': pos_array[s_idx], 'end_pos': pos_array[e_idx],
                '_s_idx': s_idx, '_e_idx': e_idx
            })
        
        if not initial_segments: continue
        df_initial = pd.DataFrame(initial_segments).sort_values('start_pos')

        # --- 3. NEW: Single-pass, ordered merging logic ---
        if df_initial.empty: continue

        final_merged_segments = []
        # Start with the first segment as the current one to merge into
        current_seg = df_initial.iloc[0].to_dict()

        for i in range(1, len(df_initial)):
            next_seg = df_initial.iloc[i].to_dict()
            gap = next_seg['start_pos'] - current_seg['end_pos']

            if gap <= merge_distance:
                # If gap is small, merge next_seg into current_seg
                current_seg['end_pos'] = max(current_seg['end_pos'], next_seg['end_pos'])
                current_seg['_e_idx'] = max(current_seg['_e_idx'], next_seg['_e_idx'])
            else:
                # If gap is too large, the current merged segment is final. Add it.
                final_merged_segments.append(current_seg)
                # The next segment becomes the new "current" segment
                current_seg = next_seg
        
        # Add the very last segment after the loop finishes
        final_merged_segments.append(current_seg)
        
        # --- 4. Final classification and stat calculation on the TRUE merged blocks ---
        final_output_rows = []
        for seg in final_merged_segments:
            s_idx, e_idx = seg['_s_idx'], seg['_e_idx']
            
            segment_states = hap_state_array[s_idx : e_idx + 1]
            n1 = np.sum(segment_states == 1)
            n2 = np.sum(segment_states == 2)
            total_snps = n1 + n2

            if total_snps < min_snps_per_segment: continue

            # Classify the final, merged segment
            label = 0
            if n1 > 0 and n2 > 0:
                minority_ratio = min(n1, n2) / total_snps
                if minority_ratio >= mosaic_minority_threshold:
                    label = 3
                else:
                    label = 1 if n1 >= n2 else 2
            elif n1 > 0: label = 1
            elif n2 > 0: label = 2
            
            if label == 0: continue

            # Calculate final probability
            prob_values = []
            indices = np.arange(s_idx, e_idx + 1)
            if label == 1:
                prob_values = prob_label1[indices[segment_states == 1]]
            elif label == 2:
                prob_values = prob_label2[indices[segment_states == 2]]
            elif label == 3:
                prob_values = np.concatenate((
                    prob_label1[indices[segment_states == 1]],
                    prob_label2[indices[segment_states == 2]]
                ))
            
            mean_prob = np.mean(prob_values) if len(prob_values) > 0 else np.nan

            final_output_rows.append({
                'chr': Chr, 'start_pos': seg['start_pos'], 'end_pos': seg['end_pos'],
                'haplotype': hap_col, 'label': label, 'snps': total_snps, 
                'prob': mean_prob, 'n_snps_label1': n1, 'n_snps_label2': n2
            })
        
        if final_output_rows:
            all_haplotypes_final_segments.append(pd.DataFrame(final_output_rows))

    if not all_haplotypes_final_segments:
        return pd.DataFrame(columns=final_output_columns)

    final_df = pd.concat(all_haplotypes_final_segments, ignore_index=True)
    return final_df[final_output_columns]
