# CHANGELOG


## v2.2.0 (2025-03-11)

### Chores

- Added examplary Dockerfile
  ([`219675f`](https://github.com/NaturalAntibody/riot_na/commit/219675f5ac80917a5604f20d66e79f29b2314d54))

### Documentation

- Corrected requirements
  ([`8741c1d`](https://github.com/NaturalAntibody/riot_na/commit/8741c1dd6c96250760b99391df68e3789bc63400))

- Update README
  ([`67f47b9`](https://github.com/NaturalAntibody/riot_na/commit/67f47b9cf8767d0683c6a794e45501301c8aca2c))

### Features

- Add get_or_create_riot_nt/aa() - cached versions of create_riot_nt/aa() with support to use in
  multiprocessing pool
  ([`e7f41a2`](https://github.com/NaturalAntibody/riot_na/commit/e7f41a278d18e63d0887bfa89a8ecd76ca67d1c1))


## v2.1.1 (2025-02-17)

### Bug Fixes

- Make int_to_str_insertion() utils method public
  ([`7b3cb9a`](https://github.com/NaturalAntibody/riot_na/commit/7b3cb9a7661389c61ad898b2433c6c3f2fa77cb7))


## v2.1.0 (2025-02-11)

### Features

- Add map insertion to number function to utils
  ([`3839e92`](https://github.com/NaturalAntibody/riot_na/commit/3839e92b514dcb09e05f5440b041ab789ef06931))


## v2.0.0 (2025-02-11)

### Features

- Add e2e tests
  ([`dbfb2d2`](https://github.com/NaturalAntibody/riot_na/commit/dbfb2d296c8f1f14645a6ffed9c549bbfff27b84))

BREAKING CHANGE: change get_region() argument from Locus to ChainType

### BREAKING CHANGES

- Change get_region() argument from Locus to ChainType


## v1.3.0 (2025-02-10)

### Build System

- Add missing maturin
  ([`445505b`](https://github.com/NaturalAntibody/riot_na/commit/445505bd1783e4b630d335e7d9388f2716f88762))

- Bump dependencies
  ([`ba5e451`](https://github.com/NaturalAntibody/riot_na/commit/ba5e45174a44d5433766a12c296f71edc43004e7))

- Migrate pyproject to poetry 2.0
  ([`8276ba6`](https://github.com/NaturalAntibody/riot_na/commit/8276ba60f4217db11efcc497486b720c84da1843))

### Features

- Add get_region() method to api
  ([`362d0b1`](https://github.com/NaturalAntibody/riot_na/commit/362d0b18ae9c758609e4c12ec4a985fb0bc3c104))


## v1.2.5 (2024-12-13)

### Continuous Integration

- Update depenendency review workflow
  ([`79f844b`](https://github.com/NaturalAntibody/riot_na/commit/79f844bce6ebea827d7e6b9ae973b8160f56f915))

### Documentation

- Add multiprocessing and Spark usage examples to README
  ([`3b0c708`](https://github.com/NaturalAntibody/riot_na/commit/3b0c70880fbef73707203192815e2e6243cce67a))


## v1.2.4 (2024-11-15)

### Bug Fixes

- Added base-ref and head-ref to dependency review action
  ([`64bfb89`](https://github.com/NaturalAntibody/riot_na/commit/64bfb89e6c2cdc06d96e9c85a69c5442511723d9))

- Importing Locus in __init__.py
  ([`4717fcd`](https://github.com/NaturalAntibody/riot_na/commit/4717fcdedd24182a40b79dd37e2178803784d6ae))

### Build System

- Bump dependencies
  ([`3f8b8e7`](https://github.com/NaturalAntibody/riot_na/commit/3f8b8e75194874f7a1c00dfc58c174e31dc298cf))

### Continuous Integration

- Add base-ref to dependency-review actions
  ([`7768fcd`](https://github.com/NaturalAntibody/riot_na/commit/7768fcda50e1661c5279ae4ae34abd880888575c))

- Added dependency review action
  ([`0007352`](https://github.com/NaturalAntibody/riot_na/commit/0007352e95dc7e549b49cc88af1ea8c9d596ea7e))


## v1.2.3 (2024-11-13)

### Bug Fixes

- Bump pyarrow version
  ([`defb606`](https://github.com/NaturalAntibody/riot_na/commit/defb6066d4867940998b6025a7e5d0680d44da69))

- Update dependencies
  ([`f724bcd`](https://github.com/NaturalAntibody/riot_na/commit/f724bcd900a659a26c9819ca2d672169c4ff4af7))

### Build System

- Pinned scikit-bio version because of danger of deprecating SSW
  ([`e3db66a`](https://github.com/NaturalAntibody/riot_na/commit/e3db66aad9d1d23a9660967b7dce1fe46b31836b))

### Refactoring

- Ignore new pylint rule
  ([`1214c86`](https://github.com/NaturalAntibody/riot_na/commit/1214c8632092627cc8f3163bba6f198cf26c9cf8))


## v1.2.2 (2024-10-24)

### Documentation

- Trigger release for docs commits
  ([`1761d9b`](https://github.com/NaturalAntibody/riot_na/commit/1761d9b11297a511ba5c90433fbf6adb6a36ec67))

### Refactoring

- Updated readme
  ([`631d630`](https://github.com/NaturalAntibody/riot_na/commit/631d630fe4ad27c6a2252f6e0e4d984f2a2b1114))


## v1.2.1 (2024-10-11)

### Bug Fixes

- Fix annotation when aa alignment failed on J gene
  ([`79d865b`](https://github.com/NaturalAntibody/riot_na/commit/79d865b601570e68e784d2f42b847fd026b942de))

### Refactoring

- Cleanup
  ([`edd5f67`](https://github.com/NaturalAntibody/riot_na/commit/edd5f6747205fbb4a66bf6db01408526b6ce52af))

- Cleanup
  ([`04bd4b8`](https://github.com/NaturalAntibody/riot_na/commit/04bd4b81a68fda7355f4e447ef5a4ebbc448de60))


## v1.2.0 (2024-10-10)

### Features

- Custom databases support
  ([`02e9dd8`](https://github.com/NaturalAntibody/riot_na/commit/02e9dd89cc64a49ed223a1e3d979f0fb9779a6ef))

- Gapped alignments
  ([`b91953a`](https://github.com/NaturalAntibody/riot_na/commit/b91953ad287d9e14b0f77461fb986a08612ee7d1))

### Refactoring

- Cleanup
  ([`572d1aa`](https://github.com/NaturalAntibody/riot_na/commit/572d1aaee93d83365ccd94bfdc1af3e511ba88a8))

- Cleanup
  ([`6cd295d`](https://github.com/NaturalAntibody/riot_na/commit/6cd295deb929c63dccfe642a50228b9ce0b883e9))


## v1.1.0 (2024-10-10)

### Features

- Added utils module with some commonly used functions
  ([`0f114b1`](https://github.com/NaturalAntibody/riot_na/commit/0f114b176b926aab7133101b0d0340abc5f42866))


## v1.0.5 (2024-09-19)

### Bug Fixes

- Unpin dependencies' fix version
  ([`8053596`](https://github.com/NaturalAntibody/riot_na/commit/8053596cacfaca9c7062b1d2180960b083ddd04c))

### Chores

- Bump scikit-bio version
  ([`058bf43`](https://github.com/NaturalAntibody/riot_na/commit/058bf436d5cb8cb2586ad1fc56a1312ca746946c))

- Regenerate poetry.lock
  ([`f280de3`](https://github.com/NaturalAntibody/riot_na/commit/f280de3dced15221c3eec21a019aa167e609da0f))

### Refactoring

- Fix linter errors
  ([`f82c6bd`](https://github.com/NaturalAntibody/riot_na/commit/f82c6bd6320c34f2c29e75071638e23b90408c16))


## v1.0.4 (2024-06-26)

### Bug Fixes

- Update README.md
  ([`8422270`](https://github.com/NaturalAntibody/riot_na/commit/84222708f7fb1c06cd74d15fb1b057d617242e2d))


## v1.0.3 (2024-06-26)

### Bug Fixes

- Update README - add links, quickstart instructions and development section
  ([`46d074e`](https://github.com/NaturalAntibody/riot_na/commit/46d074eac245ef67f3e6ce4c607b90fe94a6c907))

### Refactoring

- Added citing section
  ([`40551ee`](https://github.com/NaturalAntibody/riot_na/commit/40551eef5844eaf40587dadabaff0751441ca1e8))

- Fix workflow
  ([`a8cb445`](https://github.com/NaturalAntibody/riot_na/commit/a8cb4459f9151e98b7d203e0c94ba607d0b58feb))

- Refactored citing section
  ([`52d9c7c`](https://github.com/NaturalAntibody/riot_na/commit/52d9c7c1d6526c9458089ae27184df78867f0591))

- Stylistic change
  ([`a2f9b40`](https://github.com/NaturalAntibody/riot_na/commit/a2f9b403e3d788f368d85ac85e4f742719281798))


## v1.0.2 (2024-06-26)

### Bug Fixes

- Added short description and homepage
  ([`69985eb`](https://github.com/NaturalAntibody/riot_na/commit/69985eb3c6f26e98367636826929a725899c443b))


## v1.0.1 (2024-06-25)

### Bug Fixes

- Update README
  ([`a4d35f9`](https://github.com/NaturalAntibody/riot_na/commit/a4d35f9254c1f5edc911cd4489de326cfda5305b))


## v1.0.0 (2024-06-25)

### Features

- Init commit
  ([`8228cf2`](https://github.com/NaturalAntibody/riot_na/commit/8228cf2144a8a348c35c437121954c4797569f2d))
