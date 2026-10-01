# -*- coding: utf-8 -*-
"""final_validation.py — 投稿前最终端到端验证

一次性运行所有检查并汇总，作为可复核的完成证据。
"""
import os
import subprocess
import sys

PY_VENV = r'D:\GED_mutation\.venv\Scripts\python.exe'
PY_DSH = (r'C:\Users\黄艺林\.dsh\dsh-runtimes\dsh-primary-runtime'
          r'\dependencies\python\python.exe')
SRC = r'D:\GED_mutation\src'

CHECKS = [
    ('投稿门禁 55 项', PY_DSH, 'submission_gate.py', '总计'),
    ('稿件双向数字核验', PY_DSH, 'verify_key_numbers.py', '定向核验'),
    ('补充材料双向核验', PY_VENV, 'verify_supplementary.py', '补充材料核验'),
    ('头条数字可复现', PY_VENV, 'reproduce_headline.py', 'Recomputed'),
    ('发布包可运行性', PY_VENV, 'check_release_runnability.py', '有问题的脚本'),
    ('Methods↔代码对应', PY_DSH, 'audit_methods_code.py', 'OK'),
    ('统计方法审计', PY_DSH, 'audit_statistics.py', 'rho'),
    ('参考文献审计', PY_DSH, 'audit_references.py', '编号顺序'),
    ('取证一致性', PY_DSH, 'forensic_audit.py', '确认的内部错误'),
    ('陈旧值扫描', PY_DSH, 'final_stale_scan.py', '无陈旧值'),
    ('防御性写作', PY_DSH, 'audit_defensive.py', '需要处理'),
    ('文档健康', PY_DSH, 'audit_docx_health.py', '未引用'),
    ('图像审计', PY_DSH, 'audit_docx_images.py', 'STALE'),
]


def main():
    print('=' * 78)
    print('投稿前最终端到端验证')
    print('=' * 78)
    print()
    passed = 0
    for label, py, script, key in CHECKS:
        p = os.path.join(SRC, script)
        if not os.path.exists(p):
            print(f'  SKIP  {label:<24} (脚本不存在)')
            continue
        r = subprocess.run([py, '-u', p], capture_output=True, text=True,
                           cwd=SRC, encoding='utf-8', errors='replace')
        out = (r.stdout or '') + (r.stderr or '')
        line = next((l.strip() for l in out.split('\n') if key in l), '')
        bad = any(k in out for k in ['❌', 'FAIL', '确认的内部错误: [1-9]'])
        mark = 'OK  ' if not bad else 'FAIL'
        if not bad:
            passed += 1
        print(f'  {mark} {label:<24} {line[:58]}')

    print()
    print('=' * 78)
    print(f'端点检查：{passed}/{len(CHECKS)} 通过')
    print('=' * 78)
    return passed


if __name__ == '__main__':
    main()
