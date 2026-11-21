cd D:/Snap
PY=C:/Users/PC/miniconda3/envs/vesuvius/python.exe
$PY consistent_snap.py --windings p4/fit_A --surface p4/recto_band.zarr/2 --umbilicus p4/dl/umb_ds.json --out p4/runs/C --z0 9984 --z1 11520 --densify 2 --bbox 1180 7460 1280 7600 > p4/snap_C.log 2>&1
$PY tools/mask_flag.py p4/runs/C p4/runs/Cm 0 > p4/mask.log 2>&1
echo done > p4/snap.done
