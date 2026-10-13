"""Sphinx support for documenting CPMD input keywords.

Provides

* ``.. keyword:: NAME`` -- documents one input keyword.  Options:
  ``:section:`` (input section, e.g. ``&CPMD``), ``:options:``,
  ``:alternatives:``, ``:notes:`` (the argument summary shown next to the
  name) and ``:target:`` (link target, needed only when the same keyword
  name is documented in more than one section).
* ``:kw:`NAME``` / ``:kw:`text <TARGET>``` -- link to a keyword.
* ``.. keyword-summary:: &SECTION`` -- table of all keywords of a section.
* a "Keyword index" (``kwindex``) listing all keywords alphabetically.

The domain is registered as ``cpmd`` and is meant to be the primary domain
(``primary_domain = 'cpmd'`` in conf.py) so that the directives and roles
can be used without prefix.
"""

import re

from docutils import nodes
from docutils.parsers.rst import directives
from sphinx import addnodes
from sphinx.directives import ObjectDescription
from sphinx.domains import Domain, Index, ObjType
from sphinx.roles import XRefRole
from sphinx.util import logging
from sphinx.util.docutils import SphinxDirective
from sphinx.util.nodes import make_refnode

logger = logging.getLogger(__name__)


def make_id(target):
    return 'kw-' + re.sub(r'[^a-z0-9]+', '-', target.lower()).strip('-')


class KeywordDirective(ObjectDescription):
    has_content = True
    required_arguments = 1
    final_argument_whitespace = True
    option_spec = {
        'section': directives.unchanged_required,
        'options': directives.unchanged,
        'alternatives': directives.unchanged,
        'notes': directives.unchanged,
        'target': directives.unchanged,
        'noindex': directives.flag,
    }

    def handle_signature(self, sig, signode):
        name = sig.strip()
        signode += addnodes.desc_name(name, name)
        for opt in ('options', 'alternatives', 'notes'):
            value = self.options.get(opt)
            if value:
                signode += addnodes.desc_annotation(' ' + value, ' ' + value)
        return self.options.get('target', name).strip()

    def add_target_and_index(self, target, sig, signode):
        node_id = make_id(target)
        signode['ids'].append(node_id)
        domain = self.env.get_domain('cpmd')
        domain.add_keyword(target, sig.strip(), node_id, self.options)
        section = self.options.get('section', '')
        self.indexnode['entries'].append(
            ('single', '%s (%s keyword)' % (sig.strip(), section), node_id, '', None))

    def transform_content(self, contentnode):
        section = self.options.get('section')
        if section:
            para = nodes.paragraph()
            para += nodes.emphasis('Section: ', 'Section: ')
            para += nodes.literal(section, section)
            para['classes'].append('kw-section')
            contentnode.insert(0, para)


class KeywordIndex(Index):
    name = 'kwindex'
    localname = 'Keyword Index'
    shortname = 'keywords'

    def generate(self, docnames=None):
        content = {}
        items = sorted(self.domain.data['keywords'].items(),
                       key=lambda kv: kv[1]['name'].upper())
        for target, info in items:
            letter = info['name'][0].upper()
            content.setdefault(letter, []).append(
                (info['name'], 0, info['docname'], info['node_id'], info['section'], '', ''))
        return sorted(content.items()), True


class summary_placeholder(nodes.General, nodes.Element):
    pass


class KeywordSummary(SphinxDirective):
    required_arguments = 1
    final_argument_whitespace = True

    def run(self):
        node = summary_placeholder()
        node['section'] = self.arguments[0].strip()
        return [node]


class CPMDDomain(Domain):
    name = 'cpmd'
    label = 'CPMD'
    object_types = {'keyword': ObjType('keyword', 'kw')}
    directives = {'keyword': KeywordDirective, 'keyword-summary': KeywordSummary}
    roles = {'kw': XRefRole()}
    indices = [KeywordIndex]
    initial_data = {'keywords': {}}

    def add_keyword(self, target, name, node_id, options):
        kws = self.data['keywords']
        if target in kws and kws[target]['docname'] != self.env.docname:
            logger.warning('duplicate keyword target %r (also in %s); use :target:',
                           target, kws[target]['docname'], location=(self.env.docname, None))
        kws[target] = dict(name=name, node_id=node_id, docname=self.env.docname,
                           section=options.get('section', ''),
                           options=options.get('options', ''),
                           alternatives=options.get('alternatives', ''),
                           notes=options.get('notes', ''))

    def clear_doc(self, docname):
        for target, info in list(self.data['keywords'].items()):
            if info['docname'] == docname:
                del self.data['keywords'][target]

    def merge_domaindata(self, docnames, otherdata):
        for target, info in otherdata['keywords'].items():
            if info['docname'] in docnames:
                self.data['keywords'][target] = info

    def lookup(self, target):
        """Find a keyword by target or name; case, blanks, '_' and '-' are
        not significant, and an entry 'A, B, C' matches each of its names."""
        kws = self.data['keywords']
        target = re.sub(r'\s+', ' ', target.strip())
        if target in kws:
            return kws[target]

        def norm(x):
            return re.sub(r'[\s_\-]+', ' ', x.upper()).strip()

        key = norm(target)
        for t, info in sorted(kws.items()):
            if norm(t) == key or norm(info['name']) == key:
                return info
        for t, info in sorted(kws.items()):
            if ',' in info['name'] and key in [norm(n) for n in info['name'].split(',')]:
                return info
        return None

    def resolve_xref(self, env, fromdocname, builder, typ, target, node, contnode):
        info = self.lookup(target)
        if info is None:
            logger.warning('unknown keyword reference %r', target,
                           location=node, type='cpmd', subtype='kw')
            return None
        return make_refnode(builder, fromdocname, info['docname'], info['node_id'],
                            contnode, info['name'])

    def resolve_any_xref(self, env, fromdocname, builder, target, node, contnode):
        info = self.lookup(target)
        if info is None:
            return []
        return [('cpmd:kw', make_refnode(builder, fromdocname, info['docname'],
                                         info['node_id'], contnode, info['name']))]

    def get_objects(self):
        for target, info in self.data['keywords'].items():
            yield (target, info['name'], 'keyword', info['docname'], info['node_id'], 1)


def build_summaries(app, doctree, fromdocname):
    domain = app.env.get_domain('cpmd')
    for node in doctree.findall(summary_placeholder):
        section = node['section']
        rows = sorted((info for info in domain.data['keywords'].values()
                       if info['section'] == section),
                      key=lambda i: i['name'].upper())
        table = nodes.table(classes=['kw-summary'])
        tgroup = nodes.tgroup(cols=2)
        table += tgroup
        tgroup += nodes.colspec(colwidth=35)
        tgroup += nodes.colspec(colwidth=65)
        thead = nodes.thead()
        tgroup += thead
        hrow = nodes.row()
        for text in ('Keyword', 'Arguments'):
            entry = nodes.entry()
            entry += nodes.paragraph(text=text)
            hrow += entry
        thead += hrow
        tbody = nodes.tbody()
        tgroup += tbody
        for info in rows:
            row = nodes.row()
            entry = nodes.entry()
            para = nodes.paragraph()
            ref = make_refnode(app.builder, fromdocname, info['docname'], info['node_id'],
                               nodes.strong(info['name'], info['name']), info['name'])
            para += ref
            entry += para
            row += entry
            entry = nodes.entry()
            args = ' '.join(x for x in (info['options'], info['alternatives'], info['notes']) if x)
            entry += nodes.paragraph(text=args)
            row += entry
            tbody += row
        if not rows:
            node.replace_self(nodes.paragraph(text='(no keywords documented)'))
        else:
            node.replace_self(table)


def setup(app):
    app.add_domain(CPMDDomain)
    app.add_node(summary_placeholder)
    app.connect('doctree-resolved', build_summaries)
    return {'version': '1.0', 'parallel_read_safe': True, 'parallel_write_safe': True}
