"""Compatibility loader for Azimuth's sklearn 0.17 model pickles."""

import pickle

import numpy as np


TREE_LEAF = -1


class LeastSquaresError:
    def __setstate__(self, state):
        self.__dict__.update(state)


class MeanEstimator:
    def __setstate__(self, state):
        self.__dict__.update(state)

    def predict(self, X):
        return np.ones((X.shape[0],), dtype=np.float64) * float(self.mean)


class LegacyTree:
    def __init__(self, n_features=None, n_classes=None, n_outputs=None):
        self.n_features = n_features
        self.n_classes = n_classes
        self.n_outputs = n_outputs

    def __setstate__(self, state):
        self.max_depth = state.get("max_depth")
        self.node_count = state["node_count"]
        self.nodes = state["nodes"]
        self.values = state["values"]

    def predict(self, X):
        X = np.asarray(X)
        out = np.empty((X.shape[0],), dtype=np.float64)

        left = self.nodes["left_child"]
        right = self.nodes["right_child"]
        feature = self.nodes["feature"]
        threshold = self.nodes["threshold"]

        for row_idx in range(X.shape[0]):
            node = 0
            while left[node] != TREE_LEAF:
                if X[row_idx, feature[node]] <= threshold[node]:
                    node = left[node]
                else:
                    node = right[node]
            out[row_idx] = self.values[node].reshape(-1)[0]

        return out


class LegacyDecisionTreeRegressor:
    def __setstate__(self, state):
        self.__dict__.update(state)

    def predict(self, X):
        return self.tree_.predict(X)


class LegacyGradientBoostingRegressor:
    def __setstate__(self, state):
        self.__dict__.update(state)

    def predict(self, X):
        X = np.asarray(X)
        scores = self.init_.predict(X).astype(np.float64, copy=False)
        for stage in range(self.estimators_.shape[0]):
            tree = self.estimators_[stage, 0]
            scores += self.learning_rate * tree.predict(X)
        return scores


class _AzimuthUnpickler(pickle.Unpickler):
    _CLASS_MAP = {
        ("sklearn.ensemble.gradient_boosting", "GradientBoostingRegressor"): LegacyGradientBoostingRegressor,
        ("sklearn.ensemble.gradient_boosting", "LeastSquaresError"): LeastSquaresError,
        ("sklearn.ensemble.gradient_boosting", "MeanEstimator"): MeanEstimator,
        ("sklearn.tree.tree", "DecisionTreeRegressor"): LegacyDecisionTreeRegressor,
        ("sklearn.tree._tree", "Tree"): LegacyTree,
    }

    def find_class(self, module, name):
        mapped = self._CLASS_MAP.get((module, name))
        if mapped is not None:
            return mapped
        return super().find_class(module, name)


def load_azimuth_model(file_obj):
    return _AzimuthUnpickler(file_obj, encoding="latin1").load()
