# qm-engine-sdk

A self-contained Swift package of **63 quantum-mechanics simulation modules** with a high-performance numeric engine, extracted from the parent QM project (where the reference Python sources and the companion iOS app live). MIT licensed — see [LICENSE](LICENSE) for the "no copyright under AI" note (中文说明在 LICENSE 文件内).

- **Products**: `EngineKit` (numeric engine) + `QMModules` (simulation module library)
- **Platforms**: iOS 17+ / macOS 14+
- **Swift tools**: 6.0, Swift 6 strict concurrency
- **Dependencies**: exactly one third-party package — [`swift-numerics`](https://github.com/apple/swift-numerics) (Apple, Apache-2.0), used by the QFT module's complex channel
- **Tests**: 440+ tests across ~69 suites, all green (`swift test`)

## Installation

### Swift Package Manager (local path)

```swift
// Package.swift
dependencies: [
    .package(path: "../qm-engine-sdk"),
]
// target dependencies
.product(name: "EngineKit", package: "qm-engine-sdk"),
.product(name: "QMModules", package: "qm-engine-sdk"),
```

### Swift Package Manager (git URL — once published)

```swift
dependencies: [
    .package(url: "https://example.com/qm-engine-sdk.git", from: "1.0.0"),
]
```

### Xcode

`File → Add Package Dependencies… → Add Local…` and select this folder (or paste the git URL after publishing).

## Usage

```swift
import EngineKit
import QMModules

// Every module is a self-describing simulation unit.
let module: any SimModule = HarmonicOscillatorModule()
print(module.meta.title, module.meta.noteNumber, module.meta.tier)

// Default parameter set → compute a result.
let constants = ConstantsSet.load(.v2022)
let values = ParamValues.defaults(for: module.params)
let result = module.compute(values, constants: constants)
```

Modules advertise their metadata (`ModuleMeta`: title, category, difficulty, compute tier, notebook number), declare their parameters, and return a `SimResult` with summary values, chart-ready series, and a theory card. Compute tiers map to interaction budgets: realtime (16 ms), seconds (2 s), strategy (10 s, backed by caching + checkpointing).

## Architecture (one paragraph)

**EngineKit** is the numeric engine: `Numeric` holds `AlgebraCore` (dense/partial/tridiagonal eigensolvers over Accelerate LAPACK), FFT, and `SymEigh`; `Scheduler` provides the throttled serial compute channels that keep interactive UIs responsive; `Strategy` (StrategyCore) adds LRU result caching, convergence monitoring, and checkpoint/resume for long strategy-tier runs; `DataLayer`, `Registry`, and `Special` round out constants (CODATA), the module registry, and special functions. **QMModules** implements the 63 simulation modules (quantum foundations, wave packets, spin & qubits, solids & photons, relativity & gravity, gemology…) on top of `SimModule` protocol, each validated by physics-law unit tests against golden fixtures.

## Repository layout

```
Sources/EngineKit/     numeric engine (Numeric / Scheduler / Strategy / Registry / DataLayer / Special)
Sources/QMModules/     63 simulation modules + shared fixtures harness
Tests/EngineKitTests/  engine unit tests + fixtures
Tests/QMModulesTests/  module physics-law tests, budget tests, performance audit + fixtures
```

## License

MIT. The full text is in [LICENSE](LICENSE), which opens with a bilingual note on intent: the code is AI-generated and claims no copyright ("AI 之下无版权 / no copyright under AI") — MIT is chosen precisely because it grants maximum freedom.
