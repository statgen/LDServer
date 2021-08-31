# ldserver API

<!-- @import "[TOC]" {cmd="toc" depthFrom=2 depthTo=6 orderedList=false} -->
<!-- code_chunk_output -->

- [API specification](#api-specification)
  - [Correlation types](#correlation-types)
    - [Request](#request)
    - [Response](#response)
  - [Genome builds](#genome-builds)
    - [Request](#request-1)
    - [Response](#response-1)
  - [References](#references)
    - [Request](#request-2)
    - [Response](#response-2)
  - [Populations](#populations)
    - [Request](#request-3)
    - [Response](#response-3)
  - [Pair-wise region correlations](#pair-wise-region-correlations)
    - [Request](#request-4)
    - [Parameters](#parameters)
    - [Response](#response-4)
  - [Reference variant correlations](#reference-variant-correlations)
    - [Request](#request-5)
    - [Parameters](#parameters-1)
    - [Response](#response-5)

<!-- /code_chunk_output -->

## API specification

### Correlation types

These are the various correlation measures supported by the server. Examples are r2 and D'. 

#### Request

`GET /correlations`

#### Response

```json
{
  "data": [
    {
      "description": "",
      "label": "r",
      "name": "r",
      "type": "LD"
    },
    {
      "description": "",
      "label": "r^2",
      "name": "rsquare",
      "type": "LD"
    },
    {
      "description": "",
      "label": "covariance",
      "name": "cov",
      "type": "Covariance"
    }
  ],
  "error": null
}
```

### Genome builds

Get a list available genome builds.

#### Request

`GET /genome_builds`

#### Response

```json
{
  "data": [
    "GRCh37",
    "GRCh38"
  ],
  "error": null
}
```

### References

Get a list of all available reference panels.

#### Request

`GET /genome_builds/<build>/references`

#### Response

```json
{
  "data": [
    {
      "description": "1000 Genomes Project Phase 3",
      "name": "1000G"
    }
  ],
  "error": null
}
```

### Populations

#### Request

`GET /genome_builds/<build>/references/<reference>/populations`

#### Response

```json
{
  "data": [
    "AFR",
    "ALL",
    "AMR",
    "EAS",
    "EUR",
    "SAS"
  ],
  "error": null
}
```

### Pair-wise region correlations

#### Request

`GET /genome_builds/GRCh37/references/1000G/populations/AFR/regions`

#### Parameters

| Param | Type | Examples | Details |
| ----- | ---- | -------- | ------- | 
| chrom | string | "10", "X" |
| start | int | 1 |
| stop | int | 9 |
| correlation | string | "rsquare" | Corresponds to /correlations types
| limit | int | 1000 | Limit the number of records
| format | string | "classic" or "compact" | Classic is default. Compact is much smaller, but slightly harder to parse
| precision | int | 4 | Limits precision of LD values

#### Response

If using default format (classic):

```json
{
  "data": {
    "variant1": [
      "10:114550501_G/A",
      "10:114550501_G/A",
      "10:114550501_G/A",
      "10:114550523_G/GT",
      "10:114550523_G/GT",
      "10:114550588_CTCATAGCTGGAAGTCCCCTGGTAGGAGGA/C"
    ],
    "chromosome1": [
      "10",
      "10",
      "10",
      "10",
      "10",
      "10"
    ],
    "position1": [
      114550501,
      114550501,
      114550501,
      114550523,
      114550523,
      114550588
    ],
    "variant2": [
      "10:114550523_G/GT",
      "10:114550588_CTCATAGCTGGAAGTCCCCTGGTAGGAGGA/C",
      "10:114550600_A/G",
      "10:114550588_CTCATAGCTGGAAGTCCCCTGGTAGGAGGA/C",
      "10:114550600_A/G",
      "10:114550600_A/G"
    ],
    "chromosome2": [
      "10",
      "10",
      "10",
      "10",
      "10",
      "10"
    ],
    "position2": [
      114550523,
      114550588,
      114550600,
      114550588,
      114550600,
      114550600
    ],
    "correlation": [
      0.002376697026193142,
      0.00014225491031538695,
      0.0006637412589043379,
      0.03565419092774391,
      0.04548444598913193,
      0.0029309114906936884
    ]
  },
  "error": null,
  "next": null
}
```

If using `format=compact`:

```json
{
  "data": {
    "variants": [
      "10:114550501_G/A",
      "10:114550523_G/GT",
      "10:114550588_CTCATAGCTGGAAGTCCCCTGGTAGGAGGA/C",
      "10:114550600_A/G"
    ],
    "chromosomes": [
      "10",
      "10",
      "10",
      "10"
    ],
    "positions": [
      114550501,
      114550523,
      114550588,
      114550600
    ],
    "offsets": [
      1,
      2,
      3,
      -2147483648
    ],
    "correlations": [
      [
        0.002376697026193142,
        0.00014225491031538695,
        0.0006637412589043379
      ],
      [
        0.03565419092774391,
        0.04548444598913193
      ],
      [
        0.0029309114906936884
      ],
      []
    ]
  },
  "error": "",
  "next": ""
}
```

### Reference variant correlations

Retrieves correlations between a reference variant and all others in a defined region.

#### Request

`GET /genome_builds/GRCh37/references/1000G/populations/AFR/regions`

#### Parameters

| Param | Type | Examples | Details |
| ----- | ---- | -------- | ------- | 
| chrom | string | "10", "X" |
| start | int | 1 |
| stop | int | 9 |
| variant | string | "10:1414141_A/T" | Reference variant
| correlation | string | "rsquare" | Corresponds to /correlations types
| limit | int | 1000 | Limit the number of records
| precision | int | 4 | Limits precision of LD values

#### Response

```json
{
  "data": {
    "variant1": [
      "10:114758349_C/T",
      "10:114758349_C/T",
      "10:114758349_C/T",
      "10:114758349_C/T",
      "10:114758349_C/T",
      "10:114758349_C/T"
    ],
    "chromosome1": [
      "10",
      "10",
      "10",
      "10",
      "10",
      "10"
    ],
    "position1": [
      114758349,
      114758349,
      114758349,
      114758349,
      114758349,
      114758349
    ],
    "variant2": [
      "10:114550498_G/A",
      "10:114550501_G/A",
      "10:114550523_G/GT",
      "10:114550586_C/T",
      "10:114550588_CTCATAGCTGGAAGTCCCCTGGTAGGAGGA/C",
      "10:114550600_A/G"
    ],
    "chromosome2": [
      "10",
      "10",
      "10",
      "10",
      "10",
      "10"
    ],
    "position2": [
      114550498,
      114550501,
      114550523,
      114550586,
      114550588,
      114550600
    ],
    "correlation": [
      5.892965054954402e-05,
      2.6653198801795952e-05,
      0.0002563385060057044,
      0.0001682174624875188,
      0.0007479109917767346,
      0.0004545401025097817
    ]
  },
  "error": null,
  "next": null
}
```

If using `format=compact`:

```json
{
  "data": {
    "index_variant": "10:114758349_C/T",
    "index_chromosome": "10",
    "index_position": 114758349,
    "variants": [
      "10:114550498_G/A",
      "10:114550501_G/A",
      "10:114550523_G/GT",
      "10:114550586_C/T",
      "10:114550588_CTCATAGCTGGAAGTCCCCTGGTAGGAGGA/C",
      "10:114550600_A/G"
    ],
    "chromosomes": [
      "10",
      "10",
      "10",
      "10",
      "10",
      "10"
    ],
    "positions": [
      114550498,
      114550501,
      114550523,
      114550586,
      114550588,
      114550600
    ],
    "correlations": [
      5.892965054954402e-05,
      2.6653198801795952e-05,
      0.0002563385060057044,
      0.0001682174624875188,
      0.0007479109917767346,
      0.0004545401025097817
    ]
  },
  "error": "",
  "next": ""
}
```