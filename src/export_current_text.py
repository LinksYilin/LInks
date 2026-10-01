from docx import Document
from docx.table import Table
from docx.text.paragraph import Paragraph
P = r'D:\GED_mutation\manuscript_routeA.docx'
doc = Document(P)
lines = []
for element in doc.element.body.iterchildren():
    if element.tag.endswith('}p'):
        t = Paragraph(element,doc).text.strip()
        if t:
            lines.extend([t,''])
    elif element.tag.endswith('}tbl'):
        t = Table(element,doc)
        lines.append('[TABLE]')
        lines += [' | '.join(c.text.strip() for c in r.cells) for r in t.rows]
        lines.extend(['[/TABLE]',''])
out = r'D:\GED_mutation\manuscript_routeA_current.txt'
with open(out,'w',encoding='utf-8') as f:
    f.write('\n'.join(lines))
print(out,len(lines))
