# Quantum ESPRESSO Tutorial

An in-depth, theory + hands-on Quantum ESPRESSO tutorial, published as a
[Docusaurus](https://docusaurus.io/) website.

🌐 **Live site:** https://chaewoon11.github.io/qe-tutorial/ *(deployed via GitHub Actions)*

Every chapter pairs the **theory** with a **runnable calculation**: theory →
annotated input file → run command → output explained line by line → exercises.
Scope is core Quantum ESPRESSO (`pw.x`, DFPT, post-processing) — see the
[curriculum](docs/intro.md).

## Repository layout

| Path | What |
|---|---|
| `docs/` | the tutorial content (Markdown, rendered by Docusaurus) |
| `code/` | runnable QE inputs, reference outputs, bundled Si pseudopotential |
| `notebooks/` | Python notebooks for visualization (e.g. SCF convergence plots) |
| `static/` | images and assets served by the site |
| `docusaurus.config.mjs`, `sidebars.js` | site configuration |
| `.github/workflows/deploy.yml` | build + deploy to GitHub Pages |
| `reference/` | lab QE/EPW/HPC agent definitions, cluster guides (reference only — not part of the tutorial) |

## Develop locally

Requires Node.js ≥ 18.

```bash
npm install
npm start          # dev server with hot reload at http://localhost:3000
npm run build      # production build into build/
```

## Running the QE examples

The calculations live in `code/`. For Chapter 0:

```bash
cd code/00-first-scf
pw.x < inputs/si.scf.in > my.scf.out
```

See [`docs/setup/installation.md`](docs/setup/installation.md) for installing
Quantum ESPRESSO and pseudopotentials.

## License

Tutorial material for educational use. Quantum ESPRESSO is distributed under the
GNU GPL.
