// @ts-check

/** @type {import('@docusaurus/plugin-content-docs').SidebarsConfig} */
const sidebars = {
  tutorialSidebar: [
    'intro',
    {
      type: 'category',
      label: 'Setup',
      collapsed: false,
      items: ['setup/installation'],
    },
    {
      type: 'category',
      label: 'Tutorial',
      collapsed: false,
      items: [
        'chapters/first-scf',
        'chapters/plane-wave-dft',
        'chapters/pseudopotentials',
        'chapters/brillouin-zone',
        'chapters/relaxation',
        'chapters/bands',
        'chapters/dos',
        'chapters/charge',
        'chapters/functionals',
      ],
    },
  ],
};

export default sidebars;
