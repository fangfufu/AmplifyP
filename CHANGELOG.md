# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
