0.8.3
=====

* Set mosdepth total coverage metric to visible


0.7.0
=====

* Add tissue classifier


0.6.0
=====

* Add support for Somalier


0.5.0
=====

* Add support for mosdepth


0.4.0
=====

* Support the case that an extracted metric can have multiple types. Example: The Sequence Length reported by FastQC can either be `151` (`int`) or `10-151` (`str`).


0.3.1
=====

* Ensure all needed values are extracted from the Kraken report, set defaults when not present


0.3.0
=====

* Add support for Kraken2


0.2.1
=====

* Updated some keys and tooltips for metrics to extract


0.2.0
=====

* Add support for VerifyBamID2


0.1.0
=====

* Run Python formatter on codebase
* Added the `visible` flag to various metrics
* Add tests for bamstats
* Some refactoring to better deal with potential additions to the QualityMetric item schema
* Allow definition of additional metrics that can be calculated from extracted metrics


0.0.14
=====

* Fix typos


0.0.13
=====

* Enable RNA-SeQC


0.0.7
=====

* Minor changes to bamstats parser


0.0.3
=====

* Add FastQC


0.0.2
=====

* Add RNA-SeqQC 
* Safely cast metrics before adding them to result. If cast fails, the metric will be omitted
* Minor refactoring


0.0.1
=====

* Initial release
