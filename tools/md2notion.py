#!/usr/bin/env python3
"""Convert the derivation Markdown (GitHub/KaTeX flavour) into Notion-flavoured Markdown.
   - inline $..$ -> $`..`$ ; display $$..$$ kept as equation blocks
   - pipe tables -> <table header-row="true">
   - plain text: escape Notion specials ([ ] < > { } ^ ~ |) outside math/code
   - two-space list nesting -> tabs
   - image lines ![caption](path) kept as images; with a JSON map {path: notion source} the path is replaced
     (e.g. by the file-upload:// source returned by Notion's create-attachment)
   usage: python3 tools/md2notion.py FILE.md [IMAGES.json]"""
import re, sys, json

SPECIAL = set('[]<>{}^~|')

def esc_plain(s):
    out = []
    for ch in s:
        out.append('\\' + ch if ch in SPECIAL else ch)
    return ''.join(out)

def conv_inline(s):
    """convert a line of inline text: code spans, inline math, plain"""
    res, i, n = [], 0, len(s)
    while i < n:
        if s[i] == '`':
            j = s.find('`', i + 1)
            if j < 0: res.append(esc_plain(s[i:])); break
            res.append(s[i:j + 1]); i = j + 1
        elif s[i] == '$':
            if s.startswith('$$', i):           # inline display -> keep as inline math
                j = s.find('$$', i + 2)
                if j < 0: res.append(esc_plain(s[i:])); break
                res.append('$`' + s[i + 2:j].strip() + '`$'); i = j + 2
            else:
                j = i + 1
                while j < n and not (s[j] == '$' and s[j - 1] != '\\'): j += 1
                if j >= n: res.append(esc_plain(s[i:])); break
                res.append('$`' + s[i + 1:j] + '`$'); i = j + 1
        else:
            j = i
            while j < n and s[j] not in '`$': j += 1
            res.append(esc_plain(s[i:j])); i = j
    return ''.join(res)


def split_row(row):
    """split a pipe-table row on | that are outside $..$ math and `..` code"""
    cells, cur, inm, inc = [], [], False, False
    for ch in row:
        if ch == '`' and not inm: inc = not inc
        elif ch == '$' and not inc: inm = not inm
        if ch == '|' and not inm and not inc:
            cells.append(''.join(cur)); cur = []
        else: cur.append(ch)
    cells.append(''.join(cur))
    return [c.strip() for c in cells[1:-1]]

def indent_of(line):
    m = re.match(r'^( *)', line); return len(m.group(1))

def convert(text, images=None):
    lines = text.split('\n'); out = []; i = 0
    while i < len(lines):
        line = lines[i]
        ind = indent_of(line); tabs = '\t' * (ind // 2); body = line[ind:]
        mi = re.match(r'^!\[(.*)\]\((.*)\)\s*$', body)
        if mi:                                    # image: caption escaped, source mapped
            src = (images or {}).get(mi.group(2), mi.group(2))
            out.append(tabs + '![' + esc_plain(mi.group(1)) + '](' + src + ')'); i += 1; continue
        # display math block
        if body.startswith('$$'):
            rest = body[2:]
            if rest.strip().endswith('$$') and len(rest.strip()) > 2:      # one-line $$...$$
                eq = rest.strip()[:-2]
                out += [tabs + '$$', tabs + eq.strip(), tabs + '$$']; i += 1; continue
            eq = [rest] if rest.strip() else []
            i += 1
            while i < len(lines) and not lines[i].strip().startswith('$$'):
                eq.append(lines[i].strip()); i += 1
            i += 1
            out.append(tabs + '$$')
            out += [tabs + e for e in eq if e.strip()]
            out.append(tabs + '$$'); continue
        # table
        if body.startswith('|') and i + 1 < len(lines) and re.match(r'^\s*\|[-:| ]+\|\s*$', lines[i + 1]):
            rows = []
            while i < len(lines) and lines[i].strip().startswith('|'):
                if not re.match(r'^\s*\|[-:| ]+\|\s*$', lines[i]):
                    cells = split_row(lines[i].strip())
                    rows.append(cells)
                i += 1
            out.append('<table header-row="true">')
            for r in rows:
                out.append('\t<tr>')
                for c in r: out.append('\t\t<td>' + conv_inline(c) + '</td>')
                out.append('\t</tr>')
            out.append('</table>'); continue
        # headings / lists / quotes / plain
        m = re.match(r'^(#{1,4} |- |\d+\. |> )(.*)$', body)
        if body.strip() == '---': out.append('---')
        elif m: out.append(tabs + m.group(1) + conv_inline(m.group(2)))
        elif body.strip() == '': out.append('')
        else: out.append(tabs + conv_inline(body))
        i += 1
    return '\n'.join(out)

if __name__ == '__main__':
    src = open(sys.argv[1], encoding='utf-8').read()
    images = json.load(open(sys.argv[2])) if len(sys.argv) > 2 else None
    sys.stdout.write(convert(src, images))
