# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.8.1](https://github.com/fangfufu/AmplifyP/compare/v1.8.0...v1.8.1) (2026-07-06)


### Bug Fixes

* **ci:** use default artifact name for pages deployment ([784d4a8](https://github.com/fangfufu/AmplifyP/commit/784d4a8bcb51407dc40e80bac3697e8c60455473))

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
