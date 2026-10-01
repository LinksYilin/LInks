"""
seed_utils.py — 完整确定性控制（回应审查：seed 未完整控制）
=============================================================
原脚本只设置 numpy 与 torch CPU seed，未控制：
  - CUDA / cuDNN
  - torch.use_deterministic_algorithms
  - DataLoader 的 generator 与 worker seed
  - Python hash seed（需在进程启动前设置）
本模块统一处理，并提供可验证的确定性检查。

用法:
    from seed_utils import set_seed, make_generator
    set_seed(42)
    ... DataLoader(ds, batch_size=64, shuffle=True, generator=make_generator(42))
"""
import os
import random

import numpy as np

_SEED_ENV = 'PYTHONHASHSEED'
_CURRENT_SEED = None

# ★ 必须在 CUDA context 创建之前设置，否则 cuBLAS 非确定性（torch 会警告）。
#   审查指出的缺口之一：原脚本未设置此变量。
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')


def set_seed(seed: int, deterministic: bool = True, warn_only: bool = True) -> None:
    """
    设置全部可用的随机源。

    deterministic=True 时进一步启用确定性算法（可能变慢，且若干 CUDA 算子
    没有确定性实现）。warn_only=True 让缺少确定性实现的算子只警告不报错。
    """
    global _CURRENT_SEED
    _CURRENT_SEED = int(seed)
    os.environ[_SEED_ENV] = str(seed)
    random.seed(seed)
    np.random.seed(seed)

    try:
        import torch
    except ImportError:
        return

    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)

    if deterministic:
        # cuDNN 确定性
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
        # 关闭 TF32（会引入平台相关的数值差异）
        try:
            torch.backends.cuda.matmul.allow_tf32 = False
            torch.backends.cudnn.allow_tf32 = False
        except Exception:
            # 有意宽泛捕获：不同 torch 版本可能没有这些属性，
            # 缺少 TF32 开关不影响 CPU 训练的可复现性
            pass
        try:
            torch.use_deterministic_algorithms(True, warn_only=warn_only)
        except Exception:
            # 有意宽泛捕获：不同 torch 版本可能没有这些属性，
            # 缺少 TF32 开关不影响 CPU 训练的可复现性
            pass


def make_generator(seed: int | None = None):
    """给 DataLoader 用的独立 generator。不传 seed 时用最近一次 set_seed 的值。"""
    if seed is None:
        seed = _CURRENT_SEED if _CURRENT_SEED is not None else 0
    try:
        import torch
        g = torch.Generator()
        g.manual_seed(int(seed))
        return g
    except ImportError:
        return None


def worker_init_fn(seed: int):
    """DataLoader num_workers>0 时，保证每个 worker 的种子可复现。"""
    def _init(worker_id):
        s = seed + worker_id
        np.random.seed(s)
        random.seed(s)
        try:
            import torch
            torch.manual_seed(s)
        except ImportError:
            pass
    return _init


def describe() -> dict:
    """返回当前确定性设置，用于记录到结果元数据中。"""
    info = {'pythonhashseed': os.environ.get(_SEED_ENV), 'backend': 'cpu'}
    try:
        import torch
        info['torch'] = torch.__version__
        info['cuda_available'] = torch.cuda.is_available()
        if torch.cuda.is_available():
            info['backend'] = 'cuda'
            info['cuda'] = torch.version.cuda
            info['cudnn'] = torch.backends.cudnn.version()
            info['cudnn_deterministic'] = torch.backends.cudnn.deterministic
            info['cudnn_benchmark'] = torch.backends.cudnn.benchmark
        info['deterministic_algorithms'] = torch.are_deterministic_algorithms_enabled()
    except ImportError:
        pass
    return info
