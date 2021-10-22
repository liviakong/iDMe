for k in $(for f in `eosls $pref/Samples/AOD`; do echo $f | head -n 1 | cut -d "_" -f 2-3; done | uniq); do for t in 1 10 100 1000; do eosmkdir "$pref/Samples/AOD/${k}_ctau-${t}"; done; done

for k in $(for f in `eosls $pref/Samples/AOD/*.root`; do echo $f | head -n 1 | cut -d "_" -f 2-3; done | uniq); do for t in 1 10 100 1000; do for aod in $(eosls $pref/Samples/AOD/*${k}*ctau-${t}_*.root); do eos $cmseos mv $pref/Samples/AOD/${aod} $pref/Samples/AOD/${k}_ctau-${t}/${aod}; done; done; done
