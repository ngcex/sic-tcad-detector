from etna.core import devsim_fastsolve as _devsim_fastsolve

# devsim's UMFPACK callback without the per-solve library reload (bit-identical;
# ETNA_DEVSIM_STOCK_SOLVER=1 keeps the stock one): etna.core.devsim_fastsolve
_devsim_fastsolve.install()
