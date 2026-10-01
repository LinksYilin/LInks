# -*- coding: utf-8 -*-
"""
fix_reviewer1.py — 按审稿人 1 的确认意见修正稿件（4 项真实错误）
==================================================================
R1-M2  等价性界限是非对称的，不能声称为对称 |Δr| 边界
R1-M3  分类型接触中 formed_hydro 名义显著（p=0.0316），且未做多重比较校正
R1-M5  Cα 图并非在所有阈值下零编辑（6 Å 17 对、9/10 Å 5 对）
R1-m5  训练集为 420 个蛋白而非 88 个
"""
from docx import Document

DOC = r'D:\GED_mutation\manuscript_routeA.docx'
D = '\u0394'
M = '\u2212'
A = '\u00c5'
AL = '\u03b1'


def main():
    d = Document(DOC)
    n = 0
    repl = []

    # ---------------- R1-M2：等价性界限改为方向性表述 ----------------
    repl.append((
        'The intervals are narrower than the spread produced by varying the contact cutoff '
        '(Section 4.6, where the same quantity moves by up to 0.057 on S669 and 0.062 on ssym), so we '
        'treat them as bounding the increment rather than as a pre-specified equivalence test. On '
        'that reading they exclude structural increments larger than |' + D + 'r| \u2248 0.04 on S669 '
        'and \u2248 0.07 on ssym, and adding a contact-graph representation to a modern sequence '
        'baseline did not improve ' + D + D + 'G prediction by more than that.',
        'The intervals are asymmetric and are no narrower than the spread produced by varying the '
        'contact cutoff (Section 4.6, where the same quantity moves by up to 0.057 on S669 and 0.062 '
        'on ssym), so we report them directionally rather than as a pre-specified equivalence test. '
        'On S669 the interval excludes increments below ' + M + '0.036 and above +0.010; on ssym it '
        'excludes increments below ' + M + '0.007 and above +0.070. The negative direction on ssym is '
        'therefore only weakly constrained, and a symmetric bound cannot be quoted for either '
        'benchmark. Adding a contact-graph representation to a modern sequence baseline therefore '
        'produced no positive increment above +0.01 on S669 and none above +0.07 on ssym.'))

    repl.append((
        'which bounds structural increments below about 0.04.',
        'which on S669 excludes positive increments above +0.010 and on ssym above +0.070.'))

    # Table 1 脚注
    repl.append((
        'Read as an equivalence statement, these intervals exclude structural increments above '
        '|' + D + 'r| \u2248 0.04 on S669 and \u2248 0.07 on ssym.',
        'These intervals are asymmetric and are reported directionally: on S669 the increment is '
        'excluded below ' + M + '0.036 and above +0.010; on ssym, below ' + M + '0.007 and above '
        '+0.070. A single symmetric bound cannot be quoted for either benchmark.'))

    # ---------------- R1-M3：分类型接触的校正与准确表述 ----------------
    repl.append((
        'Stratifying by contact type did produce one significant quantity: the number of broken '
        'hydrophobic contacts correlated with ' + D + D + 'G at r = 0.135 (95% CI 0.044 to 0.201), '
        'the only contact-derived feature to exclude zero. Broken electrostatic (r = 0.047), other '
        'broken (r = 0.027), formed hydrophobic (r = 0.096), formed electrostatic (r = ' + M + '0.046) '
        'and other formed (r = 0.063) contacts were all non-significant.',
        'Stratifying by contact type produced two nominally significant quantities among six tests: '
        'the number of broken hydrophobic contacts (r = 0.135, 95% CI 0.044 to 0.201, nominal '
        'P = 0.0024) and the number of formed hydrophobic contacts (r = 0.096, nominal P = 0.032). '
        'Broken electrostatic (r = 0.047), other broken (r = 0.027), formed electrostatic '
        '(r = ' + M + '0.046) and other formed (r = 0.063) contacts were non-significant. Under a '
        'Bonferroni threshold of 0.0083 across the six tests, only the broken hydrophobic count '
        'remains significant, and we report that count with that correction in mind rather than as a '
        'feature selected from six.'))

    # ---------------- R1-M5：Cα 零编辑限定为 8 Å ----------------
    repl.append((
        AL + ' graphs recorded no contact change in any of 505 quality-controlled pairs, whereas',
        AL + ' graphs recorded no contact change in any of the 505 quality-controlled pairs at the '
        '8 ' + A + ' cutoff used throughout, whereas'))

    repl.append((
        'The C' + AL + ' graphs were identical in all 505 quality-controlled pairs, providing an '
        'internal control.',
        'The C' + AL + ' graphs were identical in all 505 quality-controlled pairs at 8 ' + A + ', '
        'providing an internal control. At 6 ' + A + ', 17 pairs (3.4%, all in protein 2PR5) differ '
        'by a single edge, and at 9 and 10 ' + A + ' five pairs (1.0%, all in 3D3B) do so; the '
        '8 ' + A + ' control is therefore exact only at the reference cutoff, and we restrict the '
        'zero-change claim to it.'))

    # ---------------- R1-m5：训练集蛋白数 ----------------
    repl.append((
        'yielding 7,905 training examples for 88 proteins.',
        'yielding 7,905 training examples drawn from 420 proteins (232 from MegaScale and 188 from '
        'ThermoMutDB).'))

    repl.append((
        'yielding 7,905 training examples for 88 proteins. S669 was used during development',
        'yielding 7,905 training examples drawn from 420 proteins (232 from MegaScale and 188 from '
        'ThermoMutDB). S669 was used during development'))

    for p in d.paragraphs:
        if not p.runs:
            continue
        full = ''.join(r.text for r in p.runs)
        new = full
        for a, b in repl:
            if a in new:
                new = new.replace(a, b)
                n += 1
        if new != full:
            p.runs[0].text = new
            for r in p.runs[1:]:
                r.text = ''

    d.save(DOC)
    print(f'已修正 {n} 处（R1-M2 / R1-M3 / R1-M5 / R1-m5）')


if __name__ == '__main__':
    main()
