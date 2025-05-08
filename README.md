# easypqp-rs

---

easypqp-rs is a Rust library for in-silico peptide library generation, with  Python bindings for integration with the python [EasyPQP](https://github.com/grosenberger/easypqp) library.

## Features

* Fast in-silico library generation using Rust

* Includes a command-line tool for batch library generation

* Python bindings for integration within the easypqp Python package

* Configurable via JSON for fine-tuning predictions, fragmentation settings, and NCE/instrument profiles

## Rust Binary CLI Example

easypqp-rs has an optional standalone command-line interface (CLI) binary for generating in-silico libraries. This can be used independently of the EasyPQP Python package if you prefer.

```bash

easypqp-insilico ./config.json
```
