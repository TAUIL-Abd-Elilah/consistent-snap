cd C:/VesuviusTemp/p4/metric
PY=C:/Users/PC/miniconda3/envs/vesuvius/python.exe
for arm in A Cm; do
  $PY get_ink_metrics.py C:/VesuviusTemp/p4/render/$arm/ink --output C:/VesuviusTemp/p4/render/$arm/metric --sequential-folds --pixel-size-um 38.4 > C:/VesuviusTemp/p4/metric_$arm.log 2>&1
done
echo done > C:/VesuviusTemp/p4/metric.done
