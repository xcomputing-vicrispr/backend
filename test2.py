import subprocess

def fold_rna(seq: str):
    result = subprocess.run(
        ["RNAfold"],
        input=seq.encode(),
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL
    )
    output = result.stdout.decode().splitlines()
    if len(output) < 2:
        return None, None
    structure = output[1].split(" ")[0]
    mfe = float(output[1].split("(")[-1].replace(")", ""))
    return structure, mfe



rev_comp = "GCGGAUUUAGCUCAGUUGG"
ss, mfe = fold_rna(rev_comp)
print("Structure:", ss)
print("MFE:", mfe)