// @ts-check
import {themes as prismThemes} from 'prism-react-renderer';
import remarkMath from 'remark-math';
import rehypeKatex from 'rehype-katex';

/** @type {import('@docusaurus/types').Config} */
const config = {
  title: 'Quantum ESPRESSO Tutorial',
  tagline: 'Theory + hands-on — from your first SCF to phonons',
  favicon: 'img/favicon.ico',

  future: {
    v4: true,
  },

  url: 'https://chaewoon11.github.io',
  baseUrl: '/qe-tutorial/',

  organizationName: 'chaewoon11',
  projectName: 'qe-tutorial',
  deploymentBranch: 'gh-pages',
  trailingSlash: false,

  onBrokenLinks: 'warn',
  onBrokenMarkdownLinks: 'warn',

  i18n: {
    defaultLocale: 'en',
    locales: ['en'],
  },

  markdown: {
    mermaid: false,
  },

  presets: [
    [
      'classic',
      /** @type {import('@docusaurus/preset-classic').Options} */
      ({
        docs: {
          routeBasePath: '/', // docs-only mode: docs served at the site root
          sidebarPath: './sidebars.js',
          editUrl: 'https://github.com/chaewoon11/qe-tutorial/tree/master/',
          remarkPlugins: [remarkMath],
          rehypePlugins: [rehypeKatex],
        },
        blog: false,
        theme: {
          customCss: './src/css/custom.css',
        },
      }),
    ],
  ],

  stylesheets: [
    {
      href: 'https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.css',
      type: 'text/css',
      integrity:
        'sha384-n8MVd4RsNIU0tAv4ct0nTaAbDJwPJzDEaqSD1odI+WdtXRGWt2kTvGFasHpSy3SV',
      crossorigin: 'anonymous',
    },
  ],

  themeConfig:
    /** @type {import('@docusaurus/preset-classic').ThemeConfig} */
    ({
      image: 'img/docusaurus-social-card.jpg',
      navbar: {
        title: 'Quantum ESPRESSO Tutorial',
        logo: {
          alt: 'QE Tutorial',
          src: 'img/logo.svg',
        },
        items: [
          {
            type: 'docSidebar',
            sidebarId: 'tutorialSidebar',
            position: 'left',
            label: 'Tutorial',
          },
          {
            href: 'https://github.com/chaewoon11/qe-tutorial',
            label: 'GitHub',
            position: 'right',
          },
        ],
      },
      footer: {
        style: 'dark',
        links: [
          {
            title: 'Tutorial',
            items: [
              {label: 'Introduction', to: '/'},
              {label: 'Setup', to: '/setup/installation'},
              {label: 'Chapter 0 — First SCF', to: '/chapters/first-scf'},
            ],
          },
          {
            title: 'Quantum ESPRESSO',
            items: [
              {label: 'Official site', href: 'https://www.quantum-espresso.org/'},
              {label: 'Documentation', href: 'https://www.quantum-espresso.org/documentation/'},
              {label: 'SSSP pseudopotentials', href: 'https://www.materialscloud.org/discover/sssp'},
            ],
          },
          {
            title: 'More',
            items: [
              {label: 'GitHub repository', href: 'https://github.com/chaewoon11/qe-tutorial'},
            ],
          },
        ],
        copyright: `Quantum ESPRESSO Tutorial. Built with Docusaurus.`,
      },
      prism: {
        theme: prismThemes.github,
        darkTheme: prismThemes.dracula,
        additionalLanguages: ['bash', 'fortran', 'python'],
      },
    }),
};

export default config;
