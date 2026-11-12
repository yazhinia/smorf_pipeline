"""Sequence hygiene gate shared by the combine step.

filterseq() is copied unchanged from classifier_dataset/rmambiguous_seq.py:32-51; the
thresholds come from assets/hygiene.yaml, which holds the hygiene block of lf_config.yaml.
"""
import re


def filterseq(cfg):
    """Return a seq -> drop-reason test; None means keep."""
    ad = cfg["hygiene"]["ambig_drop"]
    drop_x = bool(ad["drop_any_x"])
    drop_stop = bool(ad["drop_internal_stop"])
    nonstd = re.compile("[%s]" % re.escape(ad["nonstandard_residues"]))

    def keep_or_reason(seq):
        core = seq[:-1] if seq.endswith("*") else seq   # trailing stop is not internal
        if not core:
            return "empty"
        if drop_x and "X" in core:
            return "ambiguous_x"
        if drop_stop and "*" in core:
            return "internal_stop"
        if nonstd.search(core):
            return "nonstandard_residue"
        return None

    return keep_or_reason
