"""Small Biopython fallback used by the Azimuth Python 3 inference path.

The original Azimuth package depends on Biopython.  The scoring path only needs
Tm_staluc and a tiny Seq surface, so keep those pieces local for container
images that do not install Biopython.
"""

import math


_DNA_NN3 = {
    "AA": (-7.9, -22.2),
    "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7),
    "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4),
    "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0),
    "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2),
    "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9),
    "CC": (-8.0, -19.9),
}

_DNA_WEIGHTS = {
    "A": 313.21,
    "T": 304.2,
    "C": 289.18,
    "G": 329.21,
}

_COMPLEMENT = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")


def _clean(seq):
    seq = str(seq).upper().replace("U", "T")
    return "".join(base for base in seq if base in "ACGT")


class _MeltingTemp:
    @staticmethod
    def Tm_staluc(s, dnac=50, saltc=50, rna=0):
        """Return the SantaLucia nearest-neighbor Tm used by old Biopython."""
        if rna:
            raise NotImplementedError("RNA Tm_staluc fallback is not implemented")

        seq = _clean(s)
        if len(seq) < 2:
            raise ValueError("Tm_staluc requires at least two DNA bases")

        delta_h = 0.0
        delta_s = 0.0

        for base in (seq[0], seq[-1]):
            if base in "AT":
                delta_h += 2.3
                delta_s += 4.1
            else:
                delta_h += 0.1
                delta_s -= 2.8

        for i in range(len(seq) - 1):
            try:
                h, s_entropy = _DNA_NN3[seq[i : i + 2]]
            except KeyError as exc:
                raise ValueError("no thermodynamic data for neighbor %r" % seq[i : i + 2]) from exc
            delta_h += h
            delta_s += s_entropy

        delta_s += 0.368 * (len(seq) - 1) * math.log(float(saltc) * 1e-3)
        concentration = float(dnac) * 1e-9 / 4.0
        return (1000.0 * delta_h) / (delta_s + 1.987 * math.log(concentration)) - 273.15


class _SeqUtils:
    @staticmethod
    def GC(seq):
        seq = _clean(seq)
        if not seq:
            return 0.0
        return 100.0 * float(seq.count("G") + seq.count("C")) / len(seq)

    @staticmethod
    def molecular_weight(seq, seq_type="DNA"):
        seq = _clean(seq)
        if not seq:
            return 0.0
        return sum(_DNA_WEIGHTS[base] for base in seq) - 61.96


class _CompatSeq(str):
    def reverse_complement(self):
        return _CompatSeq(self.translate(_COMPLEMENT)[::-1])

    def tostring(self):
        return str(self)

    def __getitem__(self, item):
        value = super().__getitem__(item)
        if isinstance(item, slice):
            return _CompatSeq(value)
        return value


class _SeqModule:
    Seq = _CompatSeq


class _MissingModule:
    def __getattr__(self, name):
        raise ImportError("Biopython is required for Bio.%s" % name)


Tm = _MeltingTemp
SeqUtil = _SeqUtils
Seq = _SeqModule
SeqIO = _MissingModule()
Entrez = _MissingModule()

def patch_legacy_biopython(seq_util, melting_temp):
    if not hasattr(melting_temp, "Tm_staluc"):
        melting_temp.Tm_staluc = _MeltingTemp.Tm_staluc
    if not hasattr(seq_util, "GC"):
        seq_util.GC = _SeqUtils.GC
    return seq_util, melting_temp
