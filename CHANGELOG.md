# Changelog

All notable changes to this project are documented here. The format follows
[conventional commits](https://www.conventionalcommits.org); entries are
generated from commit messages by [git-cliff](https://git-cliff.org).

## [0.6.0](https://github.com/COMBINE-lab/grangers/compare/v0.5.0...v0.6.0) (2026-08-10)


### Build & CI

* Add bump_and_publish.sh and a git-cliff config ([4f203c8](https://github.com/COMBINE-lab/grangers/commit/4f203c8801fe6ce0cc35e9e3236746c932432001))

### Features

* Polars 0.53, noodles 0.115, and drop the flate2 backend pin ([a586ada](https://github.com/COMBINE-lab/grangers/commit/a586ada7288a8b4efcffab595eb245a30afed835))

### Testing

* Add parse fixtures, CI, and a real README ([9d47527](https://github.com/COMBINE-lab/grangers/commit/9d47527ce79f68b4baf2a33d7277ab54e578a4a8))

## [0.5.0](https://github.com/COMBINE-lab/grangers/compare/v0.4.0...v0.5.0) (2025-01-24)


### Bug Fixes

* GZDecoder to MultiGZDecoder ([0d00fd1](https://github.com/COMBINE-lab/grangers/commit/0d00fd1624b963c6509ea459a43d72b40bc49a4e))

### Features

* Allow reading from gzipped fasta, think more about static lifetime bound ([b235c03](https://github.com/COMBINE-lab/grangers/commit/b235c03d9430da0532960f72d56f2b9877e6d019))
* Upgrade deps ([89af5b9](https://github.com/COMBINE-lab/grangers/commit/89af5b9d6a5b7419b7d23b071086cd7e0eba0e99))

## [0.4.0](https://github.com/COMBINE-lab/grangers/compare/v0.3.1...v0.4.0) (2023-10-28)


### Bug Fixes

* Be pendantic about polars version ([2e33395](https://github.com/COMBINE-lab/grangers/commit/2e33395ee7b4a508481d47e04ae21863dcc6c32f))

### Features

* Bump noodles version ([0c93e50](https://github.com/COMBINE-lab/grangers/commit/0c93e50e46a6f8a2aa11bb2d101d8bcb9ba8652d))

## [0.3.1](https://github.com/COMBINE-lab/grangers/compare/v0.3.0...v0.3.1) (2023-07-24)


### Bug Fixes

* Merge pull request #37 from COMBINE-lab/dev ([#37](https://github.com/COMBINE-lab/grangers/issues/37))

## [0.3.0](https://github.com/COMBINE-lab/grangers/compare/v0.2.1...v0.3.0) (2023-07-24)


### Bug Fixes

* Improve parsing of fasta def in iterator ([925acd2](https://github.com/COMBINE-lab/grangers/commit/925acd2294f451c6b6131e0588cea554189ec8c1))
* Fix bugs in update_df ([142792e](https://github.com/COMBINE-lab/grangers/commit/142792e8bc02d62b79a7aca954f2fcac2bbdccbb))

### Features

* Add update_column and update_dataframe ([de48312](https://github.com/COMBINE-lab/grangers/commit/de4831270f4def5696c0c3b853ea67435757e1ad))
* Add update_column and update_df function ([8c11c10](https://github.com/COMBINE-lab/grangers/commit/8c11c109e7dd89c0cbd9738762f37e7e81b891f2))

## [0.2.1](https://github.com/COMBINE-lab/grangers/compare/v0.2.0...v0.2.1) (2023-07-10)


### Bug Fixes

* Fix cargo file ([cd17195](https://github.com/COMBINE-lab/grangers/commit/cd171953f70234d69f0f90aa1df22eda05b3364d))

## [0.2.0](https://github.com/COMBINE-lab/grangers/compare/v0.1.4...v0.2.0) (2023-07-09)


### Bug Fixes

* Format ([4cc41bd](https://github.com/COMBINE-lab/grangers/commit/4cc41bdf1e649afe121c66511373b9fa32e2282d))

### Features

* Add filterable writers ([2f4890d](https://github.com/COMBINE-lab/grangers/commit/2f4890da23d817e5e3b66b8fbd162f71ef7c7b12))

### Reformat

* Reorganize the repo to make the structure clearer ([1b1da2c](https://github.com/COMBINE-lab/grangers/commit/1b1da2c08008516cbe1d525aed4279295d083319))
* Remove lappers from Grangers as ppl can call lappers by themselves ([c9c3d39](https://github.com/COMBINE-lab/grangers/commit/c9c3d39f070884551f819557d5aba1b95be8aad8))

## [0.1.1](https://github.com/COMBINE-lab/grangers/compare/v0.1.0...v0.1.1) (2023-06-13)


### Bug Fixes

* Fix crates build ([3557af6](https://github.com/COMBINE-lab/grangers/commit/3557af6ceda2dd9695f7f85b7097e947aa6b7e8b))
* Too many keywords in crate manifest ([1708f34](https://github.com/COMBINE-lab/grangers/commit/1708f34b5e1bd71040756ffb7e1305cd42f02dba))
* Keywords in crate manifest ([8295546](https://github.com/COMBINE-lab/grangers/commit/8295546f393842ce7f3abb8485843bc95f9acdb9))
* Readme description ([43b7a79](https://github.com/COMBINE-lab/grangers/commit/43b7a7944123e6139b9770d38d7bf08cddd271bd))

## [0.1.0](https://github.com/COMBINE-lab/grangers/releases/tag/v0.1.0) (2023-06-13)


### Bug Fixes

* Fix bug in evaluation steps ([5be11d1](https://github.com/COMBINE-lab/grangers/commit/5be11d11d4768b27051fcc4ca0731856fc318710))
* Fix logging ([b681f75](https://github.com/COMBINE-lab/grangers/commit/b681f759d917b6063c537d1e79d7b20fc0ff58bb))
* Fix exon concatenation ([dc608b6](https://github.com/COMBINE-lab/grangers/commit/dc608b6831690f9559b5eef4eec768c11efed2d7))
* Fix bug for writing tx sequences ([6989ea5](https://github.com/COMBINE-lab/grangers/commit/6989ea5a41c11eaad41f7ece09eb0836eea08738))
* Made clippy happy ([a0ce59d](https://github.com/COMBINE-lab/grangers/commit/a0ce59da0a10248baf50680a30c79ff9402b7f32))
* Bug in boundary ([f36824d](https://github.com/COMBINE-lab/grangers/commit/f36824d8293784c3c7e9bf471c16a1dff97d77ef))
* Bug in boundary ([98679f0](https://github.com/COMBINE-lab/grangers/commit/98679f0172774b1b779fedb192b8c8b4e9465bd7))
* Bug in boundar ([774bbe8](https://github.com/COMBINE-lab/grangers/commit/774bbe8e814f8be5796161e3c53def91e365e4f5))
* Bug in boundary ([b01ebf1](https://github.com/COMBINE-lab/grangers/commit/b01ebf1f5289a03f35fbcaa846d42927a8952697))

### Features

* Flank function checked. ([cf54ae6](https://github.com/COMBINE-lab/grangers/commit/cf54ae6cbc5e9e34a52173bafd2855044bdffff8))
* New merge function using poloars ([a4d016e](https://github.com/COMBINE-lab/grangers/commit/a4d016eb2607d22f4d113e9093a5f9f407ca08ed))
* Add gaps option ([7058ec6](https://github.com/COMBINE-lab/grangers/commit/7058ec6d57d4c36ffa839b36878a890b81a28b77))
* Implemented introns function ([01d92fb](https://github.com/COMBINE-lab/grangers/commit/01d92fbb3f9c04275feb92b858cf8793fda36b4a))
* Testing functions and improved implementations ([4f44f8b](https://github.com/COMBINE-lab/grangers/commit/4f44f8b7fe591106765c9e37477ee66e64d5e5b7))
* Added more testing functions ([9c9bf85](https://github.com/COMBINE-lab/grangers/commit/9c9bf8504262197da8f5332b3477d605d12c6a26))
* Get_transcript_sequences function and tons of modifications to the design ([13ca1d3](https://github.com/COMBINE-lab/grangers/commit/13ca1d3cf9cf4d8bf5dbf3025d82d6fb853e32cf))
* Error checking ([13094f5](https://github.com/COMBINE-lab/grangers/commit/13094f500e219c74ad3c82f78572dac4efc4831a))

### Ffeat

* Implemented eture sequence iterator. ([7ac8332](https://github.com/COMBINE-lab/grangers/commit/7ac83320c5a36ea341e450a51b9353dd0f599138))
* Implemented eture sequence iterator. ([3395199](https://github.com/COMBINE-lab/grangers/commit/33951993e9bba08b234a793382ad333d35f14cff))


