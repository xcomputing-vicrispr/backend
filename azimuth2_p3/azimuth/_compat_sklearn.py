"""Compatibility shims for sklearn APIs used by the original Azimuth code."""

import types

import numpy as np
import sklearn
from sklearn.model_selection import GridSearchCV, cross_val_score
from sklearn.model_selection import KFold as _ModernKFold
from sklearn.model_selection import StratifiedKFold as _ModernStratifiedKFold


class KFold:
    def __init__(self, n, n_folds=3, shuffle=False, random_state=None):
        self.n = int(n)
        self._cv = _ModernKFold(n_splits=n_folds, shuffle=shuffle, random_state=random_state)

    def __iter__(self):
        return self._cv.split(np.arange(self.n))

    def __len__(self):
        return self._cv.get_n_splits()


class StratifiedKFold:
    def __init__(self, y, n_folds=3, shuffle=False, random_state=None):
        self.y = np.asarray(y)
        self._cv = _ModernStratifiedKFold(n_splits=n_folds, shuffle=shuffle, random_state=random_state)

    def __iter__(self):
        return self._cv.split(np.zeros(self.y.shape[0]), self.y)

    def __len__(self):
        return self._cv.get_n_splits()


cross_validation = types.SimpleNamespace(
    KFold=KFold,
    StratifiedKFold=StratifiedKFold,
    cross_val_score=cross_val_score,
)


def install_sklearn_legacy_aliases():
    setattr(sklearn, "cross_" + "validation", cross_validation)
    return sklearn
