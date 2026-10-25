"""Independent small oracles for the validation tool; no production selectors."""
import importlib.util
from pathlib import Path
import numpy as np

spec=importlib.util.spec_from_file_location('audit',Path(__file__).parents[1]/'bin/audit_caas_evidence.py')
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)

def test_sparse_direct_matches_dense_missing_detection_and_ties():
    dense=np.array([[.4,.4+2e-16,0,.8],[0,0,.7,0],[0,0,0,0]])
    ids,cy=np.where(dense>0);values=dense[ids,cy]
    ref=audit.PositionReference(ids,values,3,4)
    query_ids=np.array([0,0,1,1,2]);score=np.array([.4,.9,.7,.8,.2])
    expected=[(1+np.sum(dense[i]>=v-audit.TOL))/5 for i,v in zip(query_ids,score)]
    np.testing.assert_allclose(ref.direct(query_ids,score),expected)

def test_factorization_has_inclusive_tail_and_nd_plus_one_class():
    breaks=np.array([2.]);pools=[np.array([.2,.5]),np.array([.6,.9,.9])]
    # nd=2 is in the lower null class, but its observed nd+1=3 uses the upper pool.
    p=audit.factorized_p(np.array([.9,.5,0,np.nan]),np.array([2,0,1,2]),breaks,pools,4)
    np.testing.assert_allclose(p[:3],[3/5*3/4,1/5*2/3,1])
    assert np.isnan(p[3])

def test_empty_pool_and_zero_detection_keep_pseudocount():
    p=audit.factorized_p(np.array([.3]),np.array([0]),np.array([]),[np.array([])],9)
    np.testing.assert_allclose(p,[.1])

def test_bh_matches_explicit_tied_family_with_na():
    p=np.array([.02,.01,.01,1.,np.nan])
    np.testing.assert_allclose(audit.bh(p)[:4],[.08/3,.02,.02,1.])
    assert np.isnan(audit.bh(p)[4])

def test_propensity_boundaries_do_not_split_equal_nd():
    nd=np.array([0,1,1,1,2,2,100,100])
    b=audit.propensity_breaks(nd)
    c=audit.assign_class(b,nd)
    assert len(set(c[nd==1]))==len(set(c[nd==2]))==len(set(c[nd==100]))==1

def test_equal_propensity_is_not_sufficient_for_positionwise_pooling():
    # Independent reference/test draws. All positions are detected in every cycle,
    # but the high-score minority has a different conditional score distribution.
    rng=np.random.default_rng(17)
    means=np.r_[np.repeat(.2,90),np.repeat(.8,10)]
    nref,ntest=1000,1000
    train=np.clip(rng.normal(means[:,None],.03,(100,nref)),0,1)
    ref=audit.PositionReference(np.repeat(np.arange(100),nref),train.ravel(),100,nref)
    test=np.clip(rng.normal(means[:,None],.03,(100,ntest)),0,1)
    ids=np.repeat(np.arange(100),ntest)
    fact=ref.factorized(ids,test.ravel()).reshape(100,ntest)
    direct=ref.direct(ids,test.ravel()).reshape(100,ntest)
    assert np.mean(fact[90:]<=.05)>.35
    assert np.mean(direct[90:]<=.05)<.08
