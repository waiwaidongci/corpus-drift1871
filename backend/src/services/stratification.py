import json
from .binning import count_bins
from .distributions import divergence

def stratum_key(language,device):return json.dumps([language,device],ensure_ascii=False,separators=(',',':'))

def _split(rows):
 stratified=[];excluded=0
 for r in rows:
  if r['features']['language'] is None or r['features']['device'] is None:excluded+=1
  else:stratified.append(r)
 return stratified,excluded

def _groups(rows):
 out={}
 for r in rows:
  key=stratum_key(r['features']['language'],r['features']['device'])
  if key not in out:out[key]={'key':key,'language':r['features']['language'],'device':r['features']['device'],'rows':[]}
  out[key]['rows'].append(r)
 return out

def _entry(key,language,device,reference_n,candidate_n,total_reference,total_candidate,status):
 return {'key':key,'language':language,'device':device,'reference_n':reference_n,'candidate_n':candidate_n,'reference_weight':reference_n/total_reference if total_reference else 0.,'candidate_share':candidate_n/total_candidate if total_candidate else 0.,'status':status}

def build(reference,candidate,specs,features,smoothing):
 reference_valid,reference_missing=_split(reference)
 candidate_valid,candidate_missing=_split(candidate)
 reference_n=len(reference_valid);candidate_n=len(candidate_valid)
 ref_groups=_groups(reference_valid);cand_groups=_groups(candidate_valid)
 covered=sorted(set(ref_groups)&set(cand_groups),key=lambda k:(-len(ref_groups[k]['rows']),k))
 missing=sorted(set(ref_groups)-set(cand_groups),key=lambda k:(-len(ref_groups[k]['rows']),k))
 new=sorted(set(cand_groups)-set(ref_groups),key=lambda k:(-len(cand_groups[k]['rows']),k))
 entries=[]
 for key in covered+missing:
  g=ref_groups[key];c=cand_groups.get(key);entries.append(_entry(key,g['language'],g['device'],len(g['rows']),len(c['rows']) if c else 0,reference_n,candidate_n,'covered' if c else 'missing_candidate'))
 for key in new:
  g=cand_groups[key];entries.append(_entry(key,g['language'],g['device'],0,len(g['rows']),reference_n,candidate_n,'new_candidate'))
 available=bool(reference_n and candidate_n) and not missing
 new_candidate_rows=sum(e['candidate_n'] for e in entries if e['status']=='new_candidate')
 covered_candidate_rows=sum(e['candidate_n'] for e in entries if e['status']=='covered')
 missing_dimension={'reference_rows_excluded':reference_missing,'candidate_rows_excluded':candidate_missing,'note':'语言或设备为 null 的片段无法归入联合层，不参与分层加权；原始总体比较仍保留这些片段，不以零值补层。'}
 reconciliation={'reference_total_rows':len(reference),'reference_strata_sum':sum(e['reference_n'] for e in entries),'reference_missing_dimension_rows':reference_missing,'reference_balances':len(reference)==sum(e['reference_n'] for e in entries)+reference_missing,'candidate_total_rows':len(candidate),'candidate_strata_sum':sum(e['candidate_n'] for e in entries),'candidate_missing_dimension_rows':candidate_missing,'candidate_balances':len(candidate)==sum(e['candidate_n'] for e in entries)+candidate_missing,'covered_reference_weight':sum(e['reference_weight'] for e in entries if e['status']=='covered'),'covered_candidate_rows':covered_candidate_rows,'candidate_rows_in_new_strata':new_candidate_rows}
 run={'enabled':True,'dimensions':['language','device'],'bin_policy':'shared_full_reference_frozen_bins','coverage_complete':available,'weighted_available':available,'reference_stratified_n':reference_n,'candidate_stratified_n':candidate_n,'reference_strata_count':len(ref_groups),'candidate_strata_count':len(cand_groups),'covered_strata_count':len(covered),'missing_strata':[e for e in entries if e['status']=='missing_candidate'],'new_strata':[e for e in entries if e['status']=='new_candidate'],'strata':entries,'missing_dimension':missing_dimension,'reconciliation':reconciliation,'weighting_method':'post-stratification to the reference composition over observed language×device strata: reference_weight=reference_n_s/reference_stratified_n; weighted candidate bin count=sum over covered strata of candidate_count_s(b)*reference_n_s/candidate_n_s','coverage_note':None if available else ('参考没有可用的语言×设备联合层' if not reference_n else '候选缺少参考层，未计算加权分布，也不以零值补齐；缺层与新增层在分层表中单列。')}
 feature_blocks={}
 for feature in features:
  spec=specs[feature];per=[]
  for e in entries:
   rrows=ref_groups.get(e['key'],{}).get('rows',[]);crows=cand_groups.get(e['key'],{}).get('rows',[])
   rvalues=[r['features'][feature] for r in rrows];cvalues=[r['features'][feature] for r in crows]
   per.append({'key':e['key'],'language':e['language'],'device':e['device'],'reference_n':e['reference_n'],'candidate_n':e['candidate_n'],'reference_weight':e['reference_weight'],'candidate_share':e['candidate_share'],'status':e['status'],'excluded_from_weighted':e['status']!='covered','reference_missing':sum(1 for x in rvalues if x is None),'candidate_missing':sum(1 for x in cvalues if x is None),'reference_counts':count_bins(rvalues,spec),'candidate_counts':count_bins(cvalues,spec) if crows else None,'candidate_weighted_counts':None})
  block={'available':available,'reason':None if available else run['coverage_note'],'weighted_reference_counts':None,'weighted_reference_proportions':None,'weighted_candidate_counts':None,'weighted_candidate_proportions':None,'weighted_psi':None,'weighted_js':None,'weighted_review_suggested':None,'weighted_total':None,'reference_stratified_n':reference_n,'total_matches_stratified_reference':None,'tested_for_significance':False,'test_family':None,'strata':per}
  if available:
   weighted_reference_counts=[sum(row['reference_counts'][i] for row in per) for i in range(len(spec['labels']))]
   weighted_counts=[0.]*len(spec['labels'])
   for row in per:
    factor=row['reference_n']/row['candidate_n']
    for i,v in enumerate(row['candidate_counts']):weighted_counts[i]+=v*factor
    row['candidate_weighted_counts']=[v*factor for v in row['candidate_counts']]
   weighted=divergence(weighted_reference_counts,weighted_counts,smoothing)
   total=float(sum(weighted_counts))
   block.update({'weighted_reference_counts':weighted_reference_counts,'weighted_reference_proportions':weighted['reference_proportions'],'weighted_candidate_counts':weighted_counts,'weighted_candidate_proportions':weighted['candidate_proportions'],'weighted_psi':weighted['psi'],'weighted_js':weighted['js'],'weighted_review_suggested':bool(weighted['psi']>=.2 or weighted['js']>=.05),'weighted_total':total,'total_matches_stratified_reference':bool(abs(total-reference_n)<=1e-9*max(1,reference_n))})
  feature_blocks[feature]=block
 return run,feature_blocks
