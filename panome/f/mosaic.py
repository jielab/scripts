"""Module-specific literal outcome borrowing as an interpretable sensitivity."""
import numpy as np
from scipy.optimize import minimize
from borrowing import ReferenceBank


class ModuleMosaic:
    def fit(self, x, observed, ids, groups, y, w, indices, membership, names, seed, k=10):
        self.banks, self.columns, self.names = [], [], names
        self.k = min(k, len(indices))
        for module in range(len(names)):
            ix = np.flatnonzero(membership == module)
            self.columns.append(ix)
            self.banks.append(ReferenceBank(x[:, ix]/np.sqrt(len(ix)), x[:, ix], observed[:, ix],
                ids, groups, y, w, indices, seed=seed))
        self.weights = np.full(len(names), 1/len(names))
        return self

    def components(self, x, ids, groups):
        values, references = [], []
        for ix, bank in zip(self.columns, self.banks):
            risk, detail = bank.match(x[:, ix]/np.sqrt(len(ix)), ids, groups,
                                     k=self.k, temperature=1., strength=2., kernel="euclidean")
            values.append(risk)
            references.append(bank.ids[detail["jj"][:, 0]])
        return np.array(values).T, np.array(references).T

    def tune(self, values, y, w):
        weights = w/w.sum()
        def objective(beta):
            p = np.clip(values@beta, 1e-5, 1-1e-5)
            loss = np.sum(weights*(-y*np.log(p)-(1-y)*np.log1p(-p)))+.01*np.sum(beta**2)
            grad = values.T@(weights*(p-y)/(p*(1-p)))+.02*beta
            return loss, grad
        fit = minimize(objective, self.weights, jac=True, method="SLSQP",
            bounds=[(0, 1)]*len(self.weights), constraints=[{"type":"eq", "fun":lambda b:b.sum()-1,
                                                        "jac":lambda b:np.ones_like(b)}],
            options={"maxiter":1000, "ftol":1e-10})
        if not fit.success:
            raise ValueError("Module mixture failed: "+fit.message)
        self.weights = fit.x/fit.x.sum()
        return self
