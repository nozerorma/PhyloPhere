"""Dense oracle for direct block counters, including zero-based positions and ties."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest
pa = pytest.importorskip("pyarrow")
import pyarrow.parquet as pq
ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('direct_counters',ROOT/'subworkflows/SCORING/local/src/null_position_stats.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def test_direct_counts_reduce_legacy_and_new_without_conditional_pools(tmp_path):
    classes=pa.table(dict(gene=['g','g'],position=[0,1],nd=[3,1],cls=[1,2]))
    dense=np.array([[.7,.7+1e-15,.2,0],[0,.9,0,0],[0,0,0,0]])
    obs=np.array([.7,.9,.8]);genes=['g','g','h'];pos=np.array([0,1,0])
    counts=[]
    for k,cols in enumerate([[0,1],[2,3]]):
        ids,cy=np.where(dense[:,cols]>0)
        p=tmp_path/f'positions{k}.parquet'
        pq.write_table(pa.table(dict(gene=[genes[i] for i in ids],cycle=[f'b_{cols[j]+1}' for j in cy],position=pos[ids],all=dense[ids,np.array(cols)[cy]])),p)
        t=m.count_block(p,classes,np.array([1]),genes,pos,obs)
        assert t.column_names==['gene','position','k_emp']
        f=tmp_path/f'counts{k}.parquet';pq.write_table(t,f);counts.append(f)
    red=m.reduce_counts(counts,classes,genes,pos,obs,2)
    np.testing.assert_equal(red['k_emp'],(dense>=obs[:,None]-1e-12).sum(axis=1))
    np.testing.assert_equal(red['nd'],[3,1,0])
    out=tmp_path/'counters.tsv';m.write_counters(out,genes,pos,red,2)
    assert out.read_text().splitlines()[0]=='Gene\tPosition\tk_emp\tnd'
