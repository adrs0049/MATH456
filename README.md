# MATH 456 Notes

Lecture notes for MATH 456 (Mathematical Modeling) at UMass Amherst:
nonlinear dynamics for mathematical biology, from a force balance to the phase
plane. Built with [MyST](https://mystmd.org) (Jupyter Book 2).

## Build locally

```bash
git submodule update --init
pip install -r requirements.txt
make theme     # once: builds the myst-theme submodule
make start     # http://localhost:3002
make build     # static HTML in _build/html
make pdf       # PDF exports
```

## License

[CC BY-NC-SA 4.0](LICENSE)
