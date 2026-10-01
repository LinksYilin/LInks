# -*- coding: utf-8 -*-
"""fix_sup_locality.py — 补充材料：加入局域性口径修正说明"""
import os

PATHS = [r'D:\GED_mutation\补充材料_Supplementary.md',
         os.path.join(r'D:\GED_mutation', '补充材料_Supplementary.md')]

ADD = '''

---

## S10. Locality statistic: definition and sensitivity

The distance from a mutation to the contacts it loses admits more than one definition, and the reported value depends on the choice. For each mutation we measured, for each contact present in the wild-type graph and absent from the modelled mutant graph, the distance from the mutated residue to the nearest endpoint of that contact.

A contact incident on the mutated residue itself is at distance zero under this rule, and such self-contacts account for **382 of the 1,068 broken contacts** and appear in **253 of the 358** mutation pairs with at least one broken contact. Averaging over mutations without removing them gives 4.03 Å; this value is a property of the self-contact convention rather than of how far the structural change propagates.

| Definition | Broken contacts | Unchanged contacts |
|---|---|---|
| All broken contacts, per-mutation mean | 4.03 Å (95% CI 3.11–5.42) | 15.70 Å (13.95–18.18) |
| **Self-contacts removed, per-mutation mean (n = 234)** | **7.93 Å (95% CI 6.63–10.11)** | **15.93 Å (13.93–18.86)** |
| Self-contacts removed, per-contact mean | 11.93 Å | — |
| Self-contacts removed, per-mutation median | 5.88 Å | — |

The main text reports the self-contacts-removed per-mutation mean. The qualitative conclusion is unchanged under every definition: broken contacts lie markedly closer to the mutated residue than contacts that persist, and 30.1% of non-self broken contacts lie within 5 Å. The distance is not a sufficient statistic for the effect of a contact on ΔΔG, which is why the association analysis in Section 4.5 uses counts rather than distances.
'''

for p in set(PATHS):
    if not os.path.exists(p):
        continue
    t = open(p, encoding='utf-8').read()
    if 'S10. Locality statistic' in t:
        print(f'{p}: 已含 S10，跳过')
        continue
    open(p, 'w', encoding='utf-8').write(t.rstrip() + ADD)
    print(f'{p}: 已追加 S10')

t = open(PATHS[0], encoding='utf-8').read()
print()
print(f'S10 存在: {"S10. Locality statistic" in t}')
print(f'7.93 出现: {t.count("7.93")}  4.03 出现: {t.count("4.03")}')
