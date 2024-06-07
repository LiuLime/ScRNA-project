-- select top 2000 pvalue genes which p-value >=0.5 (log10p>=0.3) according to study

WITH RankedData AS (
    SELECT
        id,
        gene,
        log10p,
        abbr_id,
        ROW_NUMBER() OVER (PARTITION BY abbr_id ORDER BY ABS(log10p) DESC) AS rn
    FROM
        log10p
)
SELECT
    id,
    gene,
    log10p,
    abbr_id
FROM
    RankedData
WHERE
    rn <= 2000 and ABS(log10p)>=0.3;
