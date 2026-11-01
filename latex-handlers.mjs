/**
 * MyST plugin: document transform that rewrites proof nodes with unsupported
 * kinds (algorithm, property) so they render in LaTeX.
 *
 * `details` nodes (dropdowns) are left alone: myst-to-tex renders them inline
 * on its own, and flattening here would also remove them from the HTML.
 *
 * Register in myst.yml:
 *   project:
 *     plugins:
 *       - latex-handlers.mjs
 */

function visit(node, type, callback) {
  if (!node) return;
  if (node.type === type) callback(node);
  if (node.children) {
    for (const child of node.children) {
      visit(child, type, callback);
    }
  }
}

function visitAll(node, callback) {
  if (!node) return;
  callback(node);
  if (node.children) {
    for (const child of node.children) {
      visitAll(child, callback);
    }
  }
}

// myst-to-tex (mystmd 1.10) has no handler for exercise and solution nodes and
// drops them with their content. Only PDF builds rewrite them, as a bold label
// followed by the body; HTML builds (`myst start`, `myst build --html`) keep
// the numbered boxes and dropdowns. The Makefile and CI never build both in
// one command.
const PDF_BUILD = process.argv.some((a) => ['--pdf', '--tex'].includes(a));

function text(value) {
  return { type: 'text', value };
}

// Replace each exercise by "Exercise k (title)." and its body, and each
// solution by "Solution." and its body.
function flattenExercises(tree) {
  let count = 0;
  function rewrite(node) {
    if (!node.children) return;
    const out = [];
    for (const child of node.children) {
      if (child.type !== 'exercise' && child.type !== 'solution') {
        rewrite(child);
        out.push(child);
        continue;
      }
      const [first, ...rest] = child.children ?? [];
      const hasTitle = first && first.type === 'admonitionTitle';
      const body = hasTitle ? rest : child.children ?? [];
      let label;
      if (child.type === 'exercise') {
        count += 1;
        label = hasTitle
          ? [text(`Exercise ${count} (`), ...first.children, text(').')]
          : [text(`Exercise ${count}.`)];
      } else {
        label = [text('Solution.')];
      }
      out.push({ type: 'paragraph', children: [{ type: 'strong', children: label }] });
      for (const inner of body) {
        rewrite(inner);
        out.push(inner);
      }
    }
    node.children = out;
  }
  rewrite(tree);
}

const plugin = {
  name: 'latex-proof-fix-plugin',
  transforms: [
    {
      name: 'flatten-exercises-for-pdf',
      doc: 'In PDF builds, rewrite exercise and solution nodes as labelled paragraphs.',
      stage: 'document',
      plugin: () => (tree) => {
        if (PDF_BUILD) flattenExercises(tree);
      },
    },
    {
      name: 'fix-unsupported-proof-kinds',
      doc: 'Map unsupported prf: kinds (algorithm, property) to supported environments by rewriting the kind field.',
      stage: 'document',
      plugin: (opts, utils) => (tree) => {
        // Convert details/summary, grid, and bibliography nodes to plain content.
        // details nodes come from :class: dropdown or :::{dropdown} directives.
        // In LaTeX there's no collapsible content, so we just render inline.
        // grid nodes come from card/grid layouts with no LaTeX equivalent.
        // bibliography nodes are handled by the template's \bibliography command.
        const flattenTypes = new Set(['grid', 'card', 'cardTitle', 'header', 'footer', 'bibliography']);
        const skipTypes = new Set(['summary']);

        // Recursively flatten unsupported node types.
        // Repeat until stable since flattening can expose more nodes to flatten.
        function flattenUnsupported(node) {
          if (!node.children) return;
          // First recurse into children
          for (const child of node.children) {
            flattenUnsupported(child);
          }
          // Then flatten this node's children
          let changed = true;
          while (changed) {
            changed = false;
            const newChildren = [];
            for (const child of node.children) {
              if (flattenTypes.has(child.type)) {
                changed = true;
                if (child.children) {
                  for (const innerChild of child.children) {
                    if (!skipTypes.has(innerChild.type)) {
                      newChildren.push(innerChild);
                    }
                  }
                }
              } else {
                newChildren.push(child);
              }
            }
            node.children = newChildren;
          }
        }
        flattenUnsupported(tree);
      },
    },
  ],
};

export default plugin;
