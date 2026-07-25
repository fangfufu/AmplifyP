# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.16.0](https://github.com/fangfufu/AmplifyP/compare/v1.15.0...v1.16.0) (2026-07-25)


### Features

* **gui:** add reverse complement button for primer input ([0339cc1](https://github.com/fangfufu/AmplifyP/commit/0339cc1a1d744fc307bba7d710b1cef8e87dcc9d))
* **gui:** show binding sites when 0 amplicons ([ce7c702](https://github.com/fangfufu/AmplifyP/commit/ce7c70241bee3c20b95450500a271eee52bd403a))


### Documentation

* **changelog:** deduplicate 1.13.0 changelog entries ([3127d55](https://github.com/fangfufu/AmplifyP/commit/3127d554b437ac11b1705cb4b39a762b10e35def))
* **gui:** add Designer 1D and 2D view sections to manual ([69dce92](https://github.com/fangfufu/AmplifyP/commit/69dce92553b59924a138b575972eecd3a00227d7))
* **gui:** rewrite designer 1d view documentation ([c0da632](https://github.com/fangfufu/AmplifyP/commit/c0da6324e8820c9409276c787a9863094104beae))
* **gui:** rewrite designer 2d view documentation ([cd86f13](https://github.com/fangfufu/AmplifyP/commit/cd86f133d93ebbc2fb56458eece3c2c30c35c366))
* **gui:** rewrite GUI README and add CLI options ([f501e56](https://github.com/fangfufu/AmplifyP/commit/f501e56944cfb7864a34d9cd282f8ac0749ef0e3))
* **gui:** rewrite input view documentation ([0f26cb4](https://github.com/fangfufu/AmplifyP/commit/0f26cb48868b1edf0dedf76262979309a66b0a8c))
* **gui:** rewrite settings view user manual ([d332576](https://github.com/fangfufu/AmplifyP/commit/d33257637606bfe94fd4df6646bca5da6033968b))
* **gui:** split manual into separate view guides ([caba2b1](https://github.com/fangfufu/AmplifyP/commit/caba2b13bb823bdadcadbbe8969b1c9e3b599f62))
* **gui:** update PCR view documentation ([10e19f5](https://github.com/fangfufu/AmplifyP/commit/10e19f5cc3abc44f3e7e8567d84f18c6207f2e38))
* remove code references from GUI user manual ([7ea14b4](https://github.com/fangfufu/AmplifyP/commit/7ea14b4f5ba891084bd5ffac6994801aeb990464))
* update primer dimer view documentation ([77215b5](https://github.com/fangfufu/AmplifyP/commit/77215b50829da719ce2fe0f44fc757b67ac2f470))

## [1.15.0](https://github.com/fangfufu/AmplifyP/compare/v1.14.0...v1.15.0) (2026-07-24)


### Features

* **gui:** add save and load for designer 1d ([c4d9e5a](https://github.com/fangfufu/AmplifyP/commit/c4d9e5aee5410f12f6c1d811a3c5b9fd6de0fd96))
* **gui:** add save and load for designer 2d ([810440a](https://github.com/fangfufu/AmplifyP/commit/810440a89f9763435b5baea3bbf11b383cf70f81))
* **gui:** add screenshot export facility and window flags ([2097486](https://github.com/fangfufu/AmplifyP/commit/209748611e7d628421b4971d0737c78393f534d3))


### Bug Fixes

* **gui:** handle None parameter loading and guard min len validation ([030e98d](https://github.com/fangfufu/AmplifyP/commit/030e98dfd73a315c8fe983f3125357e99d7a8b80))
* **gui:** prevent duplicate error and support page updates in Designer 2D ([20f6c39](https://github.com/fangfufu/AmplifyP/commit/20f6c3964469c27602b15fa1681b58290cfe0f79))
* **gui:** prevent duplicate validation error in Designer 1D form ([575da49](https://github.com/fangfufu/AmplifyP/commit/575da4974d2587db1c659f15eed686269d9ae297))
* **gui:** validate all fields in Designer 2D form before returning error ([f426476](https://github.com/fangfufu/AmplifyP/commit/f4264766fd824b2ab961f89dcf9280d095143aa3))


### Documentation

* add 1D and 2D primer designer API guide ([2a2a414](https://github.com/fangfufu/AmplifyP/commit/2a2a41466a8735a3ad6797443c7067d8abfbc8c9))
* document history of Amplify ([49e895a](https://github.com/fangfufu/AmplifyP/commit/49e895a1140b922cbc2ff769b5b70e66e28561c2))
* fix GitHub alert callout rendering ([d80669d](https://github.com/fangfufu/AmplifyP/commit/d80669d646ffc0940abdcfbe896f615a5658b9f1))
* fix GitHub alert callout rendering ([a00a8ea](https://github.com/fangfufu/AmplifyP/commit/a00a8ea54f41f152447e751e5350f3cbb0a7d7e2))
* **readme:** format screenshot section and update links ([18b5429](https://github.com/fangfufu/AmplifyP/commit/18b542980e8e64cb2a9accb4ef62d1c896c30be6))
* relocate GUI manual and images to docs/gui/ ([f8b5b07](https://github.com/fangfufu/AmplifyP/commit/f8b5b07f9aef998e854165fc6e13718eaea2b289))
* revise setup and verification guides ([b4e0d03](https://github.com/fangfufu/AmplifyP/commit/b4e0d035d099cb7d4e5266e2a93ca2e3aab7bad0))
* simplify docs/README.md index ([6a8c688](https://github.com/fangfufu/AmplifyP/commit/6a8c688f8bdc26669fc556e290ece2e878c74493))
* **src:** document main.py CLI arguments ([6e28759](https://github.com/fangfufu/AmplifyP/commit/6e28759ed82b71c5e9184f832cd28619b59e6837))

## [1.14.0](https://github.com/fangfufu/AmplifyP/compare/v1.13.0...v1.14.0) (2026-07-23)


### Features

* **gui:** default 2d designer colour map to blue-orange ([332f223](https://github.com/fangfufu/AmplifyP/commit/332f223a2c7be901fa3bdb70a3c4e1f420bb40ce))
* **gui:** default 2d overlap filter to empty ([c66a59f](https://github.com/fangfufu/AmplifyP/commit/c66a59fe423d378dbdd250f3c9b0b91548adf6b8))

## [1.13.0](https://github.com/fangfufu/AmplifyP/compare/v1.12.0...v1.13.0) (2026-07-23)


### Features

* add Linux development environment setup script ([1c697d4](https://github.com/fangfufu/AmplifyP/commit/1c697d45fc19aad1c26cf7e56bf170a99c029200))
* add primer tm toggle setting and column ([fec3465](https://github.com/fangfufu/AmplifyP/commit/fec3465c644819cf62d0688f18d99581824b9513))
* add settings backup and restore functionality ([9af4a37](https://github.com/fangfufu/AmplifyP/commit/9af4a372a3ce2d8005c438b2fad736b7ba75fca1))
* **cli:** add --auto-close flag for non-interactive rendering ([d6eb7fa](https://github.com/fangfufu/AmplifyP/commit/d6eb7fa839322ce63d2cd42098307d13ed98fd0d))
* **cli:** support launching in web browser mode ([e30f1a8](https://github.com/fangfufu/AmplifyP/commit/e30f1a85826ab08c5e2f27e00f31e6b7eda3a5af))
* **dimer:** add reorder flag to dimer generation ([4c4368c](https://github.com/fangfufu/AmplifyP/commit/4c4368cb0ec780933a4606dd4112042a19834b34))
* **gui:** add 2D primer designer view ([9cb3bba](https://github.com/fangfufu/AmplifyP/commit/9cb3bbab792448b40d458a26c93994b79f6729a1))
* **gui:** add Designer 2D settings section ([a843e24](https://github.com/fangfufu/AmplifyP/commit/a843e24a5052901f101eb2f1997a7e090322963d))
* **gui:** add smooth colour temperature gradient ([309468b](https://github.com/fangfufu/AmplifyP/commit/309468b326c30d4859f30c206d0e66474a65b0f8))
* **gui:** allow PCR simulation with single primer ([5e46661](https://github.com/fangfufu/AmplifyP/commit/5e4666184b292c1a1456fb93217cf8491ee1063f))
* **gui:** configure console and file logger ([49d1184](https://github.com/fangfufu/AmplifyP/commit/49d11844c8e80549ffc546e60a47e0165f140c37))
* **gui:** configure diagnostics and logging ([f78ecdc](https://github.com/fangfufu/AmplifyP/commit/f78ecdc84e3ce7ef1a5ce638111ea63841485582))
* **gui:** group and right-align toolbar buttons ([4b03c5e](https://github.com/fangfufu/AmplifyP/commit/4b03c5e9c1e85e9235eb345bdee817802c1e17ec))
* **gui:** implement 2D results grid colour map ([2a95069](https://github.com/fangfufu/AmplifyP/commit/2a95069ec4718e81eb97924683e981a8157bf83a))
* **gui:** improve primer input layout and styling ([5d47d3e](https://github.com/fangfufu/AmplifyP/commit/5d47d3ef4dd69ca547a954bf93b85c24e22d0c14))
* **gui:** link to LICENSE in about view ([d4405a8](https://github.com/fangfufu/AmplifyP/commit/d4405a8b78221161f4d3e6687ee14530144b4447))
* **gui:** load and save only template and primers ([37a222c](https://github.com/fangfufu/AmplifyP/commit/37a222c54cb35ec0662e7b35929fcaac07a1cb8a))
* **gui:** prevent line wrapping in primer inputs ([1032a5f](https://github.com/fangfufu/AmplifyP/commit/1032a5ff7dc2d134f1ea00da2f13ba6ea08b245d))
* **gui:** update default window dimensions and 1D designer split ([2bd3505](https://github.com/fangfufu/AmplifyP/commit/2bd3505ab91466740a4d5bb118266dc7344b9582))
* **primer_designer_2d:** add 2D designer module ([ffa1f17](https://github.com/fangfufu/AmplifyP/commit/ffa1f17aebba399a08db633466f8c43aef8f42fb))
* support TSV pasting in primer textboxes ([1540576](https://github.com/fangfufu/AmplifyP/commit/154057651fd5d62a491f82214b2244cfa034f956))


### Bug Fixes

* align flet version in requirements.txt ([bb60929](https://github.com/fangfufu/AmplifyP/commit/bb60929759ec9c5d7b54e5c26a6f978d15575e6f))
* **gui:** add field-level input error indicators ([b6f5012](https://github.com/fangfufu/AmplifyP/commit/b6f5012c5421e19cbb53b335c9b90933e72b38bc))
* **gui:** fix git_dir resolution in git.py ([b46c8ed](https://github.com/fangfufu/AmplifyP/commit/b46c8edc348cd76a0b960e50d89adb83a10b543e))
* **gui:** fix pyright error in diagnostics tile ([03fdd33](https://github.com/fangfufu/AmplifyP/commit/03fdd33fc3ba7e7374a20644e6310dd500a92424))
* **gui:** fix pyright type check failure ([e2d577d](https://github.com/fangfufu/AmplifyP/commit/e2d577d48512df892584b21e7409f249da9331a3))
* **gui:** initialise logging on application startup ([9320ee2](https://github.com/fangfufu/AmplifyP/commit/9320ee22d901017bac07100f3cadeae57691ffe6))
* **gui:** prevent TSV paste from losing primer name due to UI re-extraction ([c7d594a](https://github.com/fangfufu/AmplifyP/commit/c7d594a5e8ff3139371b9e52806d622724088c3c))
* **gui:** resolve LSP type errors in designer_1d_view ([d2b1188](https://github.com/fangfufu/AmplifyP/commit/d2b11886aca0f6400560e13426b709ec263ba983))
* **gui:** set root logger to INFO to fix hang with large primer loads ([6731a35](https://github.com/fangfufu/AmplifyP/commit/6731a35a996f222bbdd78b603e0c9fc9ba945675))
* **gui:** sync dark mode state on template hot reload ([7f202d3](https://github.com/fangfufu/AmplifyP/commit/7f202d3ce3dc8bd39267a6cc01d7f2367caa54f8))
* **gui:** sync select-all checkbox cache state ([cfa80d1](https://github.com/fangfufu/AmplifyP/commit/cfa80d142589b96bf25c051506a29a10c467e7b7))
* **gui:** use equal expand for Designer 1D left/right panels ([835d880](https://github.com/fangfufu/AmplifyP/commit/835d880c446eada063e2f667abef9e46885ef0f9))
* **gui:** use Flet Colors enum for Tm text colour ([cbfec6c](https://github.com/fangfufu/AmplifyP/commit/cbfec6ca6f6890e17255d6bfc550b04c4850d76f))
* **gui:** use page.services for FilePicker in web mode ([3cdba53](https://github.com/fangfufu/AmplifyP/commit/3cdba5324212f00a2beb41e50c839b0e33abc74e))
* make bottom-left panel 2/3 of viewport height in Designer 2D view ([9674a9e](https://github.com/fangfufu/AmplifyP/commit/9674a9e75d8094082f46ddeb01f1016fb0c15904))
* make E2E dimer alignment test robust against CanvasKit OCR limits ([04395e8](https://github.com/fangfufu/AmplifyP/commit/04395e8de9d0ee576cec20afd405eeff236227f2))
* **melting:** prevent Tm division by zero ([c4f73ea](https://github.com/fangfufu/AmplifyP/commit/c4f73ea1e5d66a5b3ca878c85013ddc1c4bde953))
* remove PCR/dimer cache and fix render order ([816a604](https://github.com/fangfufu/AmplifyP/commit/816a6048230e837f34903af190a5a0bf130aafe4))
* replace deprecated Page.set_clipboard with run_task ([5d0d79c](https://github.com/fangfufu/AmplifyP/commit/5d0d79ce233c15f6fdebb390abb89228f5f238b5))
* replace float equality check with math.isclose ([e212a19](https://github.com/fangfufu/AmplifyP/commit/e212a197afa51df59e283d3ac472d1ebef8a50f8))
* replace remaining Any types with proper types in GUI modules ([7792330](https://github.com/fangfufu/AmplifyP/commit/7792330e32eadcf825a2b03517c701f50a2384b5))
* resolve 15 SonarQube quality issues ([2ad54d5](https://github.com/fangfufu/AmplifyP/commit/2ad54d5d4c0b15f81f4ca006758d951530ae2f91))
* resolve 20 SonarQube code quality issues ([6e4e7d1](https://github.com/fangfufu/AmplifyP/commit/6e4e7d1c9ef6c20723324ece2624a3101df584f6))
* resolve all float equality checks in melting.py ([af5c0ff](https://github.com/fangfufu/AmplifyP/commit/af5c0ff95f453f33d166cceee16a5fb62dd62bbc))
* resolve Codacy issues ([e4057e8](https://github.com/fangfufu/AmplifyP/commit/e4057e8ccf6800bba8b7ad14091b67e7c7517d37))
* resolve dark theme contrast and refresh lag ([0436035](https://github.com/fangfufu/AmplifyP/commit/0436035a53a37c29f535b19b97c67b2ce5432050))
* resolve LSP type errors in flet callback signatures ([2fe03fe](https://github.com/fangfufu/AmplifyP/commit/2fe03fe4c0071a4dc19f1551c60f740071c16a14))
* set left panel top/bottom split to 50/50 ([df1c69e](https://github.com/fangfufu/AmplifyP/commit/df1c69e30be92bf45b57ed8b4076235e57841744))
* **setup:** address PR review feedback ([b8d9314](https://github.com/fangfufu/AmplifyP/commit/b8d931425d80ddfc1d96a6115a20a65ad68ff154))
* show Mean Overlap instead of Min Overlap in 2D result card ([c3f2ca3](https://github.com/fangfufu/AmplifyP/commit/c3f2ca364bc5f4afd97b40695b5011383ae6c80b))
* stretch initial 2D results grid prompt to fill available space ([2511a74](https://github.com/fangfufu/AmplifyP/commit/2511a7497f0c20359b583280eb1370a0f8efb006))
* strip newline on Enter key press in input fields ([5e9e567](https://github.com/fangfufu/AmplifyP/commit/5e9e567704ffae2288f9d8e3e02847f55f73a2ca))
* **tests:** resolve relative path in desktop e2e test ([15cac0e](https://github.com/fangfufu/AmplifyP/commit/15cac0ef9aa1f473c1ac3ad6c9b0bac5a7099150))
* update error banner colors on theme change ([ce099e2](https://github.com/fangfufu/AmplifyP/commit/ce099e266fe82edcbb7c293c4287d19ff2a6e02c))
* update Flet event handler type annotations ([b0327be](https://github.com/fangfufu/AmplifyP/commit/b0327be5f55a584b9af810c7d3eaef48b831fdff))
* use theme-aware colors for error banner ([9cb16f4](https://github.com/fangfufu/AmplifyP/commit/9cb16f4809c622fce577544ca11ef9aebfec7e2f))


### Performance Improvements

* **dimer:** Optimize PrimerDimerGenerator double loop ([c78ba74](https://github.com/fangfufu/AmplifyP/commit/c78ba74bdd05bd213f550f9fa20233ab2e4e50f9))
* **gui:** optimise save state loading ([7a3e3b1](https://github.com/fangfufu/AmplifyP/commit/7a3e3b164307c0d42a8203dc13c474f9f1d1a430))
* Optimize `Repliconf.search` inner loop lookups ([36ca432](https://github.com/fangfufu/AmplifyP/commit/36ca432eee034131428120a2a718e4b65a1b29f0))
* optimize BasePairWeightsTbl initialization to reduce list allocations ([28d9cd3](https://github.com/fangfufu/AmplifyP/commit/28d9cd3268cce67d949f4933c3f75c448900cd52))


### Documentation

* add binary download info to README ([f5d3046](https://github.com/fangfufu/AmplifyP/commit/f5d30460b81b92fa4916a25c4c21775cdbe3f63f))
* link Amplify4 and Roboto Mono in README.md ([7023a14](https://github.com/fangfufu/AmplifyP/commit/7023a146f30e03b731d00687af1013b890c75f8c))
* **readme:** prevent wrapping release-please link ([2ce1d38](https://github.com/fangfufu/AmplifyP/commit/2ce1d38b6c7a887f569d4c502594f80873aa586d))
* remove deprecated repository badges from README ([03c2f96](https://github.com/fangfufu/AmplifyP/commit/03c2f9658ea9c4481cb17426b32a69d24e88e114))
* rewrite and expand all documentation ([a6dfcdc](https://github.com/fangfufu/AmplifyP/commit/a6dfcdca8f3dce60a776d5832648bd0ad4b70d4f))
* simplify Amplify4 attribution wording ([7de6a66](https://github.com/fangfufu/AmplifyP/commit/7de6a66ab5b76354cd0d6506909de8549a59a208))


### Miscellaneous Chores

* force release version 1.13.0 ([1e64582](https://github.com/fangfufu/AmplifyP/commit/1e64582d89fb91f899631a14f363cc1f489d654a))

## [1.12.0](https://github.com/fangfufu/AmplifyP/compare/v1.11.0...v1.12.0) (2026-07-22)


### Features

* **gui:** add 1D primer designer view ([58c3497](https://github.com/fangfufu/AmplifyP/commit/58c34977fcaee307258fc5cc82d7b595983dc4a4))
* **primer_designer:** add quality_score and best_score to PrimerDesigner1D ([6c3faf2](https://github.com/fangfufu/AmplifyP/commit/6c3faf263954fe98cec53a665cbbcc951b7379d8))


### Bug Fixes

* remove unnecessary isinstance checks in PrimerDesigner1D ([8a3f380](https://github.com/fangfufu/AmplifyP/commit/8a3f3808ad45d8cecb99036f9512cd8aa5a38a68))


### Documentation

* update and add docstrings to improve code documentation ([1a63a0a](https://github.com/fangfufu/AmplifyP/commit/1a63a0a87bebbdf2fb6445dc540f5df855675ac3))
* update outdated class docstrings across codebase ([916456b](https://github.com/fangfufu/AmplifyP/commit/916456b17661051dcce33caccc4d17ab6cdb19a8))

## [1.11.0](https://github.com/fangfufu/AmplifyP/compare/v1.10.1...v1.11.0) (2026-07-22)


### Features

* **primer:** add PrimerDesigner1D module ([f5fad68](https://github.com/fangfufu/AmplifyP/commit/f5fad68c7d0fea690050a31e2fe7879776a20526))


### Bug Fixes

* address PR review feedback for PR [#302](https://github.com/fangfufu/AmplifyP/issues/302) ([4f86061](https://github.com/fangfufu/AmplifyP/commit/4f860619b893c7715a04295c1ef2a46f30bcd2bc))
* images in gui_guide.md ([e9e67cf](https://github.com/fangfufu/AmplifyP/commit/e9e67cf5c2c47d947ab4b618aee3a6128802d797))
* resolve pyright errors in primer designer ([8f77e07](https://github.com/fangfufu/AmplifyP/commit/8f77e072af2fa44867956aefde5ae19849c4f6ae))
* resolve SonarQube issues in PR [#302](https://github.com/fangfufu/AmplifyP/issues/302) ([c8c5f04](https://github.com/fangfufu/AmplifyP/commit/c8c5f041ec374dfc84f099ead12fd4df7dc4c830))


### Documentation

* rewrite README and update documentation screenshots ([492c935](https://github.com/fangfufu/AmplifyP/commit/492c93589ce4230301b18fe5e9106aa392051be3))

## [1.10.1](https://github.com/fangfufu/AmplifyP/compare/v1.10.0...v1.10.1) (2026-07-21)


### Bug Fixes

* **gui:** guard page resolution against detached control ([f6579d3](https://github.com/fangfufu/AmplifyP/commit/f6579d3ee86b300be1072b1ac34a4b692cfc2374))
* **gui:** preserve cursor on arrow navigation ([d7ac914](https://github.com/fangfufu/AmplifyP/commit/d7ac914f7581931cc509e454a3387f2badf177b3))

## [1.10.0](https://github.com/fangfufu/AmplifyP/compare/v1.9.0...v1.10.0) (2026-07-20)


### Features

* add left/right arrow navigation between primer fields ([30890a1](https://github.com/fangfufu/AmplifyP/commit/30890a1f89ce60ea765cbbce23f471368a214668))
* **gui:** add Auto wrap option to selector ([6f24ea3](https://github.com/fangfufu/AmplifyP/commit/6f24ea3cab565a638e9c4aba71c3b00a9619c462))
* **gui:** add auto-activate setting and primer list ([f4945fb](https://github.com/fangfufu/AmplifyP/commit/f4945fbb4fe86ca125e58bd925ee3e5ddcc7e6fc))
* **gui:** default bases per line to Auto ([eb90da3](https://github.com/fangfufu/AmplifyP/commit/eb90da395225e1c438f13d107e42d1b158d2d44c))
* **gui:** display self-dimer in primer info box ([5fb42f8](https://github.com/fangfufu/AmplifyP/commit/5fb42f86d0a9b49b3a1947b465b61eb02f81381b))
* **gui:** enable keyboard navigation in primer input panel ([2cfe1c2](https://github.com/fangfufu/AmplifyP/commit/2cfe1c284521a742dcacbe8bcb88890830737c66))
* **gui:** keep primer info box visible in fixed height mode ([7ebec64](https://github.com/fangfufu/AmplifyP/commit/7ebec643834a909b87e0ac219ffb1e8f906ec58f))
* **gui:** reset primer row scroll on blur ([639af1b](https://github.com/fangfufu/AmplifyP/commit/639af1b0b61e23197cdbd31b48f5cdac8f1cdd35))
* **gui:** restrict wrap length options to 10-100 ([fab55ff](https://github.com/fangfufu/AmplifyP/commit/fab55ff941e36b3dcfe5bbae92e727987dba6c1f))
* **gui:** separate header buttons with vertical divider ([3c6361a](https://github.com/fangfufu/AmplifyP/commit/3c6361ae1f495f240cac82424d498c4990d74edd))
* **gui:** show template scrollbar always at top ([eb06a8b](https://github.com/fangfufu/AmplifyP/commit/eb06a8b0c1221c61c39db67253c0989a985d977d))
* **gui:** support repositioning primer info panel ([1539690](https://github.com/fangfufu/AmplifyP/commit/15396908ff16b1713cb76ad5f6f1bd9e651af2d3))
* **gui:** update dimer view primer name alignment ([7e9ce57](https://github.com/fangfufu/AmplifyP/commit/7e9ce57eac1f555f21028648ed6c6d98c0ea28be))
* improve cursor positioning and Tab navigation in primer fields ([c08b333](https://github.com/fangfufu/AmplifyP/commit/c08b3339fe3db530774d6f782a22e8e27637fec8))


### Bug Fixes

* add isinstance check for ft.TextField before accessing .value ([5c5eab5](https://github.com/fangfufu/AmplifyP/commit/5c5eab55499ec91837e9880581f8aec4f473ba49))
* codacy issue ([d7480cb](https://github.com/fangfufu/AmplifyP/commit/d7480cbb66e8bd23ea6094c65276e98a85a7ade9))
* **gui:** batch checkbox updates via page.update() for semantics ([4c096a3](https://github.com/fangfufu/AmplifyP/commit/4c096a34aa9d6e8a9f2895dc2a4f22fdfb1f6a30))
* **gui:** calculate template panel width dynamically ([c42e426](https://github.com/fangfufu/AmplifyP/commit/c42e42648a80b9a220477c0d94b8bf870336ab73))
* **gui:** disable autowrap in template input ([63b75c0](https://github.com/fangfufu/AmplifyP/commit/63b75c040d54c0f97847b88345bb66378804c3e9))
* **gui:** ensure template sequence formatting and gutter reflow on all changes ([b54372d](https://github.com/fangfufu/AmplifyP/commit/b54372db8a39e083477ffa96d8369bf2790e231b))
* **gui:** expand template text field to full width ([22939ce](https://github.com/fangfufu/AmplifyP/commit/22939ceba79394036d652b0a549b8e860ce43e7d))
* **gui:** fix primer auto-activation race in e2e ([144aff2](https://github.com/fangfufu/AmplifyP/commit/144aff2b6d82b177fa0a925aa836714d3faa2af3))
* **gui:** increase circular checkbox container width to 115px ([33d6db0](https://github.com/fangfufu/AmplifyP/commit/33d6db0c8e45dce7789a8a386b4f9b72abef1cf9))
* **gui:** increase safety margin and char width for template wrap ([15c2039](https://github.com/fangfufu/AmplifyP/commit/15c2039a683b1683dcf56d5684ee2a12cbe4bd12))
* **gui:** increase template wrap safety margin to 100px ([f6ed585](https://github.com/fangfufu/AmplifyP/commit/f6ed585961ffe9220fc01a7cb6f2e06324a05745))
* **gui:** update template gutter base markers immediately on paste ([86e80a4](https://github.com/fangfufu/AmplifyP/commit/86e80a4f5d60e7f573518508a4f0d15b5a2c2a9b))
* Various fixes ([f314dc0](https://github.com/fangfufu/AmplifyP/commit/f314dc037c2524b3ecd48cf261958486c6c8c50b))

## [1.9.0](https://github.com/fangfufu/AmplifyP/compare/v1.8.1...v1.9.0) (2026-07-16)


### Features

* **deps:** add ignore list and comment support ([41be8ab](https://github.com/fangfufu/AmplifyP/commit/41be8ab6326ef43894f85e13c3d7fd5dafe08a31))


### Bug Fixes

* **deps:** address pull request review comments ([bf35bc5](https://github.com/fangfufu/AmplifyP/commit/bf35bc57c77def68282a76c6e7b44f7bd90ffe4f))


### Performance Improvements

* **deps:** cache package version requests to PyPI ([827242d](https://github.com/fangfufu/AmplifyP/commit/827242de609f47563b322eb06d91ab480e81502a))

## [1.8.1](https://github.com/fangfufu/AmplifyP/compare/v1.8.0...v1.8.1) (2026-07-16)


### Bug Fixes

* **web:** fix pyodide startup crash and pin flet ([7f5850e](https://github.com/fangfufu/AmplifyP/commit/7f5850e17aebf6a32ad4e7134421819d3f796467))

## [1.8.0](https://github.com/fangfufu/AmplifyP/compare/v1.7.0...v1.8.0) (2026-07-06)


### Features

* **gui:** add fixed width option for template sequence ([87fd508](https://github.com/fangfufu/AmplifyP/commit/87fd508aad28b69a4f0baeb2fbb958731acd9ea9))


### Bug Fixes

* **gui:** enforce bases per line range of 10 to 10000 ([3dc58a8](https://github.com/fangfufu/AmplifyP/commit/3dc58a862196204a640d82fd60b74a9d5bace21d))

## [1.7.0](https://github.com/fangfufu/AmplifyP/compare/v1.6.0...v1.7.0) (2026-07-06)


### Features

* **gui:** add clear all button with confirmation ([5474911](https://github.com/fangfufu/AmplifyP/commit/5474911dfc51e8860c2035a5bc2c51bdca754051))
* **gui:** add template sequence casing buttons ([4d368cc](https://github.com/fangfufu/AmplifyP/commit/4d368cc31dcac04f669a81403b81263b8007289f))


### Bug Fixes

* **gui:** address additional review comments ([1ecd69e](https://github.com/fangfufu/AmplifyP/commit/1ecd69e88b1e9e24a6e4b7e4355140b09ea64cff))
* **gui:** address PR review comments ([c332e4a](https://github.com/fangfufu/AmplifyP/commit/c332e4a586492c338ed587bc97b05b0dd3c76228))
* **gui:** fix clear dialogue pyright errors ([69d25bc](https://github.com/fangfufu/AmplifyP/commit/69d25bc713d4d89026801daa92fc1e7d6307e5c7))
* **gui:** keep clear dialogue in overlay ([dd2d2c1](https://github.com/fangfufu/AmplifyP/commit/dd2d2c1451dc5c4b1638a71aa2eeb06f055c81bd))
* **gui:** remove duplicate change_selection_case ([2e23f71](https://github.com/fangfufu/AmplifyP/commit/2e23f717d74c468cda57f12683f2e9f9aed92bbc))

## [1.6.0](https://github.com/fangfufu/AmplifyP/compare/v1.5.0...v1.6.0) (2026-07-06)


### Features

* **gui:** add status bar to template input ([04d683e](https://github.com/fangfufu/AmplifyP/commit/04d683e78680869246236fb9ed55921b14042501))
* **gui:** show total bases on template blur ([adcbaa2](https://github.com/fangfufu/AmplifyP/commit/adcbaa2f4dbc5787849922e5074db472174266ab))


### Bug Fixes

* **gui:** allow system theme without snapping ([bee7d04](https://github.com/fangfufu/AmplifyP/commit/bee7d0420b020eebcdaa80994b66698fe50924c0))
* **gui:** guard against zero font size ([731e0b2](https://github.com/fangfufu/AmplifyP/commit/731e0b2068bd64fd1b982957cfa015657d5413de))
* **gui:** handle dark_mode coercion before bool check ([e324e0e](https://github.com/fangfufu/AmplifyP/commit/e324e0ed7494a1b63f6dd4f3754cdeace1bf4422))
* **gui:** resolve Pyright errors in TemplateInput ([2cbbc84](https://github.com/fangfufu/AmplifyP/commit/2cbbc8418e251647620db79fd2c2199b18dff6e3))
* ignore readonly textareas in PRIMER_INPUT_SEL ([a97f0bf](https://github.com/fangfufu/AmplifyP/commit/a97f0bfa072ef9b7a20efa83a1b5862688b28b4e))
* **input_view:** catch RuntimeError when accessing page in timer_callback ([aa95a5a](https://github.com/fangfufu/AmplifyP/commit/aa95a5a788283a6b5e15ee9ce18a33246e96002e))


### Performance Improvements

* **gui:** optimize drag and gutter calculation ([290a68a](https://github.com/fangfufu/AmplifyP/commit/290a68a174ecd34be489aac28a1a3ba116f74231))
* **gui:** optimize layout and sequence updates ([53b9974](https://github.com/fangfufu/AmplifyP/commit/53b99745005f585c70f93af7f745ecd39ac4b03f))

## [1.5.0](https://github.com/fangfufu/AmplifyP/compare/v1.4.0...v1.5.0) (2026-07-05)


### Features

* **gui:** add auto-reload for template and primers ([ade6d36](https://github.com/fangfufu/AmplifyP/commit/ade6d363054e17cc9a7b9def5ee50ac0bb4d0e94))
* **installer:** add shortcut options ([5358521](https://github.com/fangfufu/AmplifyP/commit/53585218594980ebe3c83704bb93e3b9d151e48d))
* **installer:** use native Start Menu selection ([32eb9bd](https://github.com/fangfufu/AmplifyP/commit/32eb9bd71189807d1864d035096ce4e13c28694d))


### Bug Fixes

* **gui:** fix low contrast selected rows ([db72e3d](https://github.com/fangfufu/AmplifyP/commit/db72e3de50b675a98066b3766e3602947814ad05))


### Performance Improvements

* **gui:** optimize Flet UI updates and reduce lag ([c375234](https://github.com/fangfufu/AmplifyP/commit/c375234c7251db15c71756691d193723d5452980))
* **gui:** split on_change and stop_editing ([46d2a47](https://github.com/fangfufu/AmplifyP/commit/46d2a475e0aad704fbe8ad4f59b14534b0eaf68a))


### Documentation

* add README index page for docs subfolder ([9e8f91a](https://github.com/fangfufu/AmplifyP/commit/9e8f91ad9ae9ab45a8454920fb8537a9fce68e08))
* add Windows development setup guide ([658a1a3](https://github.com/fangfufu/AmplifyP/commit/658a1a3886b413f082df048904120f0ec8a42b7a))

## [1.4.0](https://github.com/fangfufu/AmplifyP/compare/v1.3.2...v1.4.0) (2026-07-04)


### Features

* add configurable background version checking ([32832fc](https://github.com/fangfufu/AmplifyP/commit/32832fc3f522b08384299238d49e870417fe6911))
* add Inno Setup installer build for Windows ([79dc79d](https://github.com/fangfufu/AmplifyP/commit/79dc79d275bf93e8441c358534abfab5f6205cee))


### Bug Fixes

* **build:** correct src path in gen_git_sha ([2424843](https://github.com/fangfufu/AmplifyP/commit/2424843059438373aad7bb64ff50a75fb15242ae))
* **ci:** fallback to pyproject.toml version ([aed7515](https://github.com/fangfufu/AmplifyP/commit/aed7515fc5b604f43df2adb8b93d40ed29b28c2e))
* **ci:** make pages artifact name unique ([c3020cb](https://github.com/fangfufu/AmplifyP/commit/c3020cbea96d0544953d5998ffbb9437469bd02f))
* **ci:** remove backticks from installer version ([a87c405](https://github.com/fangfufu/AmplifyP/commit/a87c40580eb15300bcf578bc4b622bd6539318c4))
* **ci:** use icon.ico for installer icon ([58fb80b](https://github.com/fangfufu/AmplifyP/commit/58fb80b5d52dad00da633e61e004728cf693f80e))
* **gui:** fix mixed-precision version comparison ([5546151](https://github.com/fangfufu/AmplifyP/commit/554615129c99d4aad45beebc0cb111be4afd6ee7))
* use [[ instead of [ for shell conditionals (SonarQube shelldre:S7688) ([96ec98a](https://github.com/fangfufu/AmplifyP/commit/96ec98a1de9baff581c3529dbb9a93655fa4504e))

## [1.3.2](https://github.com/fangfufu/AmplifyP/compare/v1.3.0...v1.3.2) (2026-07-03)


### Bug Fixes

* toggle Show Primer Temperature in List without restart ([5486903](https://github.com/fangfufu/AmplifyP/commit/5486903e7e6e3831fb68fbe7da98b08a77647755))


### Miscellaneous Chores

* force release ([bbdd844](https://github.com/fangfufu/AmplifyP/commit/bbdd844e304b03d3ac16bbba0c787cb58c641e00))

## [1.3.0](https://github.com/fangfufu/AmplifyP/compare/v1.2.0...v1.3.0) (2026-07-03)


### Features

* add double-click range selection for primer rows ([bec49e5](https://github.com/fangfufu/AmplifyP/commit/bec49e54ecdc3651eecf4c1a893b3c311ed91724))

## [1.2.0](https://github.com/fangfufu/AmplifyP/compare/v1.1.0...v1.2.0) (2026-07-01)


### Features

* add settings to ignore inactive primers in duplicate checks ([8cc286b](https://github.com/fangfufu/AmplifyP/commit/8cc286bc11cd63d26cbdf301d40a414aa48d0cf5))
* **appearance:** increase dropdown width from 350 to 500 ([3e11cc8](https://github.com/fangfufu/AmplifyP/commit/3e11cc81085cc5c49ba36ef7c94c7104020b4088))
* **diagnostics:** use BorderedCheckbox and increase dropdown width to 500 ([d2afcdc](https://github.com/fangfufu/AmplifyP/commit/d2afcdc7c3180f130b86ad895d9a643241d7e50c))
* **gui:** add configurable width to ScoreTable ([bf9519a](https://github.com/fangfufu/AmplifyP/commit/bf9519ab547d11e32e6cfb4bda9d1dfebfb3fc2d))
* **score:** increase dropdown width to 500 and set table width to 810 ([cdc4354](https://github.com/fangfufu/AmplifyP/commit/cdc43543cfc0fd4eb9b1fc3f4aa3f4005d552bfe))
* **tm:** increase dropdown width from 450 to 500 ([18d0dac](https://github.com/fangfufu/AmplifyP/commit/18d0dac5fc7ebf27016f75719c7ad84f78f9df0d))
* wire up ignore_inactive settings across input layer ([a57332c](https://github.com/fangfufu/AmplifyP/commit/a57332c264eb9d1f8db90451984eb0b71e4695d3))


### Bug Fixes

* **gui:** remove stale duplicate primer check ([48b2d73](https://github.com/fangfufu/AmplifyP/commit/48b2d737e55c5342425c65ce75026f6877ab9a7d))
* narrow None checks in primer_row and update typos hook ([610b3d2](https://github.com/fangfufu/AmplifyP/commit/610b3d2e47a569586254d329f5ae258fdf566d3b))


### Performance Improvements

* batch-update rows instead of individual row.update() ([07fa4ff](https://github.com/fangfufu/AmplifyP/commit/07fa4ff9a4695cd771c7ffa3a1910c5fead05858))

## [1.1.0](https://github.com/fangfufu/AmplifyP/compare/v1.0.1...v1.1.0) (2026-07-01)


### Features

* add live drag-and-drop row reordering ([57c818f](https://github.com/fangfufu/AmplifyP/commit/57c818f15949299ab08da84250d6e7d8e98a1b50))
* implement multi-select and decouple selection ([f669c47](https://github.com/fangfufu/AmplifyP/commit/f669c478414de178dddd84c4bb1fdd32aee3c2e9))
* improve paste behavior ([fd51ece](https://github.com/fangfufu/AmplifyP/commit/fd51eceb09a0fcfe6e89d6dedea3d3e3d87f355a))


### Bug Fixes

* resolve out-of-bounds crash when copying primers ([9bd5f9d](https://github.com/fangfufu/AmplifyP/commit/9bd5f9d455d080b74a61001f0408ce04b7024fec))

## [1.0.1](https://github.com/fangfufu/AmplifyP/compare/v1.0.0...v1.0.1) (2026-06-30)

* **chore:** Various refactorings.

## [1.0.0](https://github.com/fangfufu/AmplifyP/commit/19d347b9a3efce223385f54777efd97694919022) - 2026-06-30

Initial release
