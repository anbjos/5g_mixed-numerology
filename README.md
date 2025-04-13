Absolutely — here’s a unified and polished introduction that combines your motivation for mixed numerologies## 5g_mixed-numerology

This repository contains a subsystem-level model addressing a central challenge in [5G](https://en.wikipedia.org/wiki/5G) systems: how to schedule and process [mixed numerologies](https://www.rfwireless-world.com/articles/5g/5g-nr-numerology-terminology).

While 4G was primarily designed to **connect people**, 5G expands that vision to **connect everything**—from autonomous vehicles navigating busy roads to battery-powered water meters in deep basements. These devices have vastly different physical characteristics and communication needs. High-mobility use cases require low-latency and high-reliability links, whereas low-power sensors may only need to transmit a few bits infrequently. Supporting such diverse requirements calls for a more flexible physical layer.

To achieve this, 5G introduces the concept of **multiple numerologies**—sets of OFDM parameters (like subcarrier spacing, symbol duration, and cyclic prefix) that can coexist within the same system. However, enabling multiple numerologies side by side introduces a non-trivial challenge: **Inter-Numerology Interference (INI)**.

### Inter-Numerology Interference (INI)

Each numerology defines its own OFDM symbol duration. Since OFDM relies on the assumption that each subcarrier is a periodic signal over the symbol time, this becomes problematic when one numerology observes or processes another using its own (mismatched) symbol length. The signal is no longer sampled over an integer number of periods, which leads to **discontinuities** in the time domain and **spectral leakage** in the frequency domain.

This violation of periodicity breaks the orthogonality between subcarriers and causes **Inter-Numerology Interference**—even if numerologies are strictly separated in time or frequency. INI is thus not just a result of overlap, but a structural issue caused by differing symbol lengths and timing assumptions across numerologies.

Addressing INI requires more than just guard bands or simple filtering. It demands careful scheduling, timing-aware resource allocation, and potentially architectural separation of numerology domains. This repository explores one such scheduling concept at a high level of abstraction. The model is designed to demonstrate feasibility, with an emphasis on **modularity and flexibility**—making it a useful tool for **architectural exploration** and **system-level optimization**, rather than low-level implementation.

### Concept of the Solution

The proposed solution addresses Inter-Numerology Interference (INI) through a modular signal processing pipeline that treats each numerology independently, while enabling efficient and interference-aware mixing when appropriate. The core ideas are:

1. **Band Limitation via Filtering**  
   Each numerology is constrained to a well-defined frequency band using a low-pass filter. This filter suppresses out-of-band energy caused by discontinuities between symbols, minimizing spectral leakage and reducing INI.

2. **Per-Numerology Filtering at Low Sample Rates**  
   To reduce computational overhead, filtering is applied individually to each numerology at its **native (low) sample rate**—i.e., the rate matched to its symbol duration and subcarrier spacing. This enables efficient isolation and pre-processing without unnecessary upsampling.

3. **Interpolation for Mixing and Resource Sharing**  
   After filtering, each numerology is **interpolated to a common higher sample rate**. Although the sampling resolution increases, each signal’s bandwidth remains narrow. This creates an opportunity to **combine multiple numerologies** at the same sample rate, provided their total bandwidth fits within the available spectrum and they are **spectrally adjacent**.

4. **Local Mixing Decisions Based on Sample Rate Compatibility**  
   The decision to mix numerologies is made **locally** at each stage, using only data from numerologies that already share the same sample rate. If mixing isn't feasible at a given stage (due to overlap or bandwidth constraints), numerologies are **interpolated to the next higher sample rate**, and the mixing decision is reevaluated.

5. **Final Baseband Stage at Output Sample Rate**  
   At the highest sample rate—used for baseband processing—all remaining numerologies are **fully mixed into a single composite signal**. This is possible because, by this point, all signals have been filtered, resampled, and allocated in frequency such that **no spectral overlap occurs**, assuming that frequency planning by the network has been done correctly. This final stage prepares the unified baseband signal for transmission or further downstream processing.

