# piMIMIC 0.3.0

# piMIMIC 0.2.0

# piMIMIC 0.3.0 (2026-08-20)

- Removed dependency on the `scripty` package.
- Now uses `semTools::indProd` directly to create product indicators.
- Fixed a bug in the `$DIF.Global` output that duplicated results for each item. Now each item appears only once.
- Improved internal documentation and added comments in English to facilitate understanding of the process.

# piMIMIC 0.2.0 (2026-07-20)

- Renamed from piRFA to piMIMIC.
- Implements the `piMIMIClrt()` function for piMIMIC DIF analysis using likelihood ratio tests.
- Presentation improvements.

# piMIMIC 0.1.0 (2024-05-10)

- Initial package version.
- Implements the `piMIMIC()` function for DIF analysis using the Product Indicators (PI) approach and the Score test.
- Includes Oort's adjustment to control Type I error rate.
