This contains an experimental partial version of the Taxonium backend implemented in Rust.

Behaviour:
- without string interning it is significantly faster than the JS backend (but uses lots of memory)
- [with string interning](https://github.com/theosanderson/taxrust/pull/1) it is still a bit faster than JS, but not extremely different. It also maybe still uses more memory than the JS version - we should figure that out.

With WASM we would have a 4 GB memory limit for the web worker, which is problematic.

For now we will not proceed with this.
