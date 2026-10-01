# ✅ DOI 已完成

> **已完成。** 概念 DOI：10.5281/zenodo.23086956（始终指向最新版）
> 版本 DOI（v1.0.1）：10.5281/zenodo.23086957
> 归档页：https://zenodo.org/records/23086957
> 仓库：https://github.com/LinksYilin/LInks
> 原始操作步骤保留在下方，仅供追溯。

---

# 如何拿到 DOI —— 只剩你点 3 下

> **我已经完成的**：
> - ✅ 289 个文件已推送到 https://github.com/LinksYilin/LInks
> - ✅ 已创建 Release **v1.0.0**：https://github.com/LinksYilin/LInks/releases/tag/v1.0.0
>
> **我为什么不能代做最后一步**：DOI 由 Zenodo 铸造，而 Zenodo 要求**用你本人的账号登录一次**（OAuth 授权）。这一步需要浏览器点击"授权"，任何程序都无法代你完成——它代表你本人同意把代码归档到你的名下。
>
> **你只需点 3 下**，然后告诉我，剩下的我全做。

---

## 你的 3 下（约 2 分钟）

### 第 1 下：登录 Zenodo

打开 → **https://zenodo.org**

点右上角 **Sign up**（或 **Log in**）

选择 **Sign up with GitHub** / **Log in with GitHub**

> 如果出现"Authorize Zenodo"授权页，点 **Authorize**。这是 Zenodo 要读取你的仓库列表，属正常。

---

### 第 2 下：打开集成设置

登录后直接打开 → **https://zenodo.org/account/settings/github/**

在列表里找到 **LInks**

把它左边的开关**拨到 ON**

> 如果列表里没有 LInks，点 **Sync now** 刷新一下。

---

### 第 3 下：告诉我

回来说一句「**好了**」即可。

**我会立刻**：删除并重建 v1.0.0 Release → Zenodo 的 webhook 触发 → 自动生成 DOI → 我把 DOI 填进论文的 4 个位置。

---

## 之后会发生什么

1. Zenodo 把你的代码归档，生成形如 `10.5281/zenodo.1234567` 的 DOI
2. Zenodo 会给你的 GitHub 仓库自动加一个 badge
3. 我查询 Zenodo API 拿到 DOI，填入：
   - 摘要 Availability
   - Data & Code Availability
   - 补充材料 S9
   - 标题页 / CITATION.cff
4. 重跑全部门禁与核验
5. **投稿就绪**

---

## 如果你想自己看结果

归档成功后，这两个地址能查到：

- 你的 Zenodo 记录：https://zenodo.org/account/settings/github/
- 搜索你的仓库：https://zenodo.org/search?q=LinksYilin

---

## 备选方案（如果 GitHub 集成有问题）

Zenodo 也支持**手动上传**：

1. 打开 → https://zenodo.org/uploads/new
2. **Upload** 上传 `D:\GED_mutation\release` 打包的 zip（我可以帮你打包）
3. 填 Title / Authors / Description（我可以把内容准备好给你复制）
4. 点 **Publish** → 生成 DOI

需要我打包 zip 并生成可直接粘贴的元数据，说一声即可。

---

## 当前状态

| 项 | 状态 |
|---|---|
| GitHub 仓库 | ✅ 289 文件已推送 |
| Release v1.0.0 | ✅ 已创建 |
| Zenodo 集成 | ⏳ **等你开启** |
| DOI | ⏳ 开启后自动生成 |
| 论文 | ✅ 22 项缺陷已修，75 项数字核验通过 |
| 作者信息 | ⏳ 等你提供 |
