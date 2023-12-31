WITH pop AS
(
SELECT 
       AGG_STU_ENRL_MVW.PIDM,
       AGG_STU_ENRL_MVW.TRANSFER_CHRT_TRM_CODE
FROM EDW.AGG_STU_ENRL_MVW AGG_STU_ENRL_MVW
WHERE (AGG_STU_ENRL_MVW.TERM_CODE >= '201405')
      AND (AGG_STU_ENRL_MVW.LEVEL_NAME_EDW_DEFINITION = 'Undergraduate')
      AND (AGG_STU_ENRL_MVW.HDCNT_OFFCL_INCL_IND = 'Y')
      AND (AGG_STU_ENRL_MVW.STU_TYPE_CODE IN ('C', 'F')
      AND (AGG_STU_ENRL_MVW.DECEASED_IND != 'Y'))
)

SELECT DISTINCT
       pop.PIDM,
       FIRST_TERM,
       CURRENT_PROGRAM,
       RESIDENCY,
       CASE WHEN M1113TERM <= FIRST_TERM AND M1113TERM IS NOT NULL THEN 1 ELSE 0 END AS "MATH1113",
       CASE WHEN M1551TERM <= FIRST_TERM AND M1551TERM IS NOT NULL THEN 1 ELSE 0 END AS "MATH1551",
       CASE WHEN M1552TERM <= FIRST_TERM AND M1552TERM IS NOT NULL THEN 1 ELSE 0 END AS "MATH1552",
       CASE WHEN M1553TERM <= FIRST_TERM AND M1553TERM IS NOT NULL THEN 1 ELSE 0 END AS "MATH1553",
       CASE WHEN M1554TERM <= FIRST_TERM AND M1554TERM IS NOT NULL THEN 1 ELSE 0 END AS "MATH1554",
       CASE WHEN M2550TERM <= FIRST_TERM AND M2550TERM IS NOT NULL THEN 1 ELSE 0 END AS "MATH2550",
       CASE WHEN M2551TERM <= FIRST_TERM AND M2551TERM IS NOT NULL THEN 1 ELSE 0 END AS "MATH2551",
       CASE WHEN CS1301TERM <= FIRST_TERM AND CS1301TERM IS NOT NULL THEN 1 ELSE 0 END AS "CS1301",
       CASE WHEN CS1331TERM <= FIRST_TERM AND CS1331TERM IS NOT NULL THEN 1 ELSE 0 END AS "CS1331",
       CASE WHEN CS1371TERM <= FIRST_TERM AND CS1371TERM IS NOT NULL THEN 1 ELSE 0 END AS "CS1371",
       CASE WHEN PHYS2211TERM <= FIRST_TERM AND PHYS2211TERM IS NOT NULL THEN 1 ELSE 0 END AS "PHYS2211",
       CASE WHEN PHYS2212TERM <= FIRST_TERM AND PHYS2212TERM IS NOT NULL THEN 1 ELSE 0 END AS "PHYS2212",
       CASE WHEN START_PROGRAM <> CURRENT_PROGRAM THEN 1 ELSE 0 END AS "Changed_Major",
       CASE WHEN GRADTERM IS NOT NULL THEN GRADTERM ELSE Null END AS "Grad_Term",
       GPA,
       HRS_EARNED
FROM pop

--gets application data (not currently using)
LEFT OUTER JOIN
(
  SELECT 
    a.PIDM AS APIDM,
    a.APPLICATION_TERM_CODE AS APPTERM
  FROM EDW.AGG_STU_ADMIT_MVW a 
  WHERE a.APPLICANT_TYPE_CODE = 'F'
  AND a.ENROLLED = '1'
  AND a.LATEST_DCSN_CODE LIKE 'A%'
) ON APIDM = pop.PIDM

--gets FIRST TERM as First-year
LEFT OUTER JOIN
(SELECT
    ENRL.PIDM AS FPIDM,
    ENRL.PRI_PRGM_CODE AS START_PROGRAM,
    MIN(ENRL.TERM_CODE) AS FIRST_TERM		
FROM		
	EDW.AGG_STU_ENRL_MVW ENRL
WHERE 
     ENRL.STU_TYPE_CODE = 'F'
     AND ENRL.SNAPSHOT_IND = 'Y' --ADDED
GROUP BY ENRL.PIDM, ENRL.PRI_PRGM_CODE
) ON FPIDM = pop.PIDM AND FIRST_TERM IS NOT Null

--gets GPA data as undergrad
LEFT OUTER JOIN
(SELECT g2.PIDM AS GPIDM,
g2.TRM_CODE AS GTERM,
g2.LVL_GPA AS GPA,
g2.LVL_GPA_HRS_EARNED AS HRS_EARNED
FROM (
       SELECT g.*, row_number() OVER (PARTITION BY g.PIDM ORDER BY g.TRM_CODE DESC) AS rn
       FROM EDW.AGG_STU_GPA_TERM_MVW g WHERE g.LVL_CODE = 'US' AND g.GPA_TYPE_CODE = 'I'
     ) g2
WHERE g2.rn = 1
) ON GPIDM = pop.PIDM

--gets MAX TERM info
LEFT OUTER JOIN
(SELECT
    ENRL.PIDM AS MTPIDM,
    ENRL.TERM_CODE AS MTTERM,
    ENRL.PRI_DEGR_CODE AS CURRENT_PROGRAM,
    ENRL.RESIDENCY_NAME AS RESIDENCY
FROM (
       SELECT e.*, row_number() OVER (PARTITION BY e.PIDM ORDER BY e.TERM_CODE DESC) AS rn
       FROM EDW.AGG_STU_ENRL_MVW e WHERE e.LEVEL_CODE = 'US' AND e.PRI_COLLEGE_CODE != 'R'
     ) ENRL
WHERE ENRL.rn = 1
) ON MTPIDM = pop.PIDM

--checks for MATH1113 pre-enrollment credit
LEFT OUTER JOIN
(
SELECT DISTINCT 
s1.SHRTRCE_PIDM AS M1113PIDM,
s1.SHRTRCE_TERM_CODE_EFF AS M1113TERM
FROM EDWSTG.STG_BAN_SHRTRCE_TBL s1
INNER JOIN EDWSTG.STG_BAN_SHRTRCR_TBL s2
ON s1.SHRTRCE_PIDM = s2.SHRTRCR_PIDM AND s1.SHRTRCE_TRCR_SEQ_NO = s2.SHRTRCR_SEQ_NO
INNER JOIN EDWSTG.STG_BAN_SHRTRIT_TBL s3
ON s3.SHRTRIT_PIDM = s2.SHRTRCR_PIDM AND s2.SHRTRCR_TRIT_SEQ_NO = s3.SHRTRIT_SEQ_NO
WHERE s1.SHRTRCE_SUBJ_CODE = 'MATH'
AND s1.SHRTRCE_CRSE_NUMB = '1113'
AND s1.SHRTRCE_CREDIT_HOURS > 0
) ON M1113PIDM = pop.PIDM

--checks for MATH1551 pre-enrollment credit
LEFT OUTER JOIN
(
SELECT DISTINCT 
s1.SHRTRCE_PIDM AS M1551PIDM,
s1.SHRTRCE_TERM_CODE_EFF AS M1551TERM
FROM EDWSTG.STG_BAN_SHRTRCE_TBL s1
INNER JOIN EDWSTG.STG_BAN_SHRTRCR_TBL s2
ON s1.SHRTRCE_PIDM = s2.SHRTRCR_PIDM AND s1.SHRTRCE_TRCR_SEQ_NO = s2.SHRTRCR_SEQ_NO
INNER JOIN EDWSTG.STG_BAN_SHRTRIT_TBL s3
ON s3.SHRTRIT_PIDM = s2.SHRTRCR_PIDM AND s2.SHRTRCR_TRIT_SEQ_NO = s3.SHRTRIT_SEQ_NO
WHERE s1.SHRTRCE_SUBJ_CODE = 'MATH'
AND s1.SHRTRCE_CRSE_NUMB = '1551'
AND s1.SHRTRCE_CREDIT_HOURS > 0
) ON M1551PIDM = pop.PIDM

--checks for MATH1552 pre-enrollment credit
LEFT OUTER JOIN
(
SELECT DISTINCT 
s1.SHRTRCE_PIDM AS M1552PIDM,
s1.SHRTRCE_TERM_CODE_EFF AS M1552TERM
FROM EDWSTG.STG_BAN_SHRTRCE_TBL s1
INNER JOIN EDWSTG.STG_BAN_SHRTRCR_TBL s2
ON s1.SHRTRCE_PIDM = s2.SHRTRCR_PIDM AND s1.SHRTRCE_TRCR_SEQ_NO = s2.SHRTRCR_SEQ_NO
INNER JOIN EDWSTG.STG_BAN_SHRTRIT_TBL s3
ON s3.SHRTRIT_PIDM = s2.SHRTRCR_PIDM AND s2.SHRTRCR_TRIT_SEQ_NO = s3.SHRTRIT_SEQ_NO
WHERE s1.SHRTRCE_SUBJ_CODE = 'MATH'
AND s1.SHRTRCE_CRSE_NUMB = '1552'
AND s1.SHRTRCE_CREDIT_HOURS > 0
) ON M1552PIDM = pop.PIDM

--checks for MATH1553 pre-enrollment credit
LEFT OUTER JOIN
(
SELECT DISTINCT 
s1.SHRTRCE_PIDM AS M1553PIDM,
s1.SHRTRCE_TERM_CODE_EFF AS M1553TERM
FROM EDWSTG.STG_BAN_SHRTRCE_TBL s1
INNER JOIN EDWSTG.STG_BAN_SHRTRCR_TBL s2
ON s1.SHRTRCE_PIDM = s2.SHRTRCR_PIDM AND s1.SHRTRCE_TRCR_SEQ_NO = s2.SHRTRCR_SEQ_NO
INNER JOIN EDWSTG.STG_BAN_SHRTRIT_TBL s3
ON s3.SHRTRIT_PIDM = s2.SHRTRCR_PIDM AND s2.SHRTRCR_TRIT_SEQ_NO = s3.SHRTRIT_SEQ_NO
WHERE s1.SHRTRCE_SUBJ_CODE = 'MATH'
AND s1.SHRTRCE_CRSE_NUMB = '1553'
AND s1.SHRTRCE_CREDIT_HOURS > 0
) ON M1553PIDM = pop.PIDM

--checks for MATH1554 pre-enrollment credit
LEFT OUTER JOIN
(
SELECT DISTINCT 
s1.SHRTRCE_PIDM AS M1554PIDM,
s1.SHRTRCE_TERM_CODE_EFF AS M1554TERM
FROM EDWSTG.STG_BAN_SHRTRCE_TBL s1
INNER JOIN EDWSTG.STG_BAN_SHRTRCR_TBL s2
ON s1.SHRTRCE_PIDM = s2.SHRTRCR_PIDM AND s1.SHRTRCE_TRCR_SEQ_NO = s2.SHRTRCR_SEQ_NO
INNER JOIN EDWSTG.STG_BAN_SHRTRIT_TBL s3
ON s3.SHRTRIT_PIDM = s2.SHRTRCR_PIDM AND s2.SHRTRCR_TRIT_SEQ_NO = s3.SHRTRIT_SEQ_NO
WHERE s1.SHRTRCE_SUBJ_CODE = 'MATH'
AND s1.SHRTRCE_CRSE_NUMB = '1554'
AND s1.SHRTRCE_CREDIT_HOURS > 0
) ON M1554PIDM = pop.PIDM

--checks for MATH2550 pre-enrollment credit
LEFT OUTER JOIN
(
SELECT DISTINCT 
s1.SHRTRCE_PIDM AS M2550PIDM,
s1.SHRTRCE_TERM_CODE_EFF AS M2550TERM
FROM EDWSTG.STG_BAN_SHRTRCE_TBL s1
INNER JOIN EDWSTG.STG_BAN_SHRTRCR_TBL s2
ON s1.SHRTRCE_PIDM = s2.SHRTRCR_PIDM AND s1.SHRTRCE_TRCR_SEQ_NO = s2.SHRTRCR_SEQ_NO
INNER JOIN EDWSTG.STG_BAN_SHRTRIT_TBL s3
ON s3.SHRTRIT_PIDM = s2.SHRTRCR_PIDM AND s2.SHRTRCR_TRIT_SEQ_NO = s3.SHRTRIT_SEQ_NO
WHERE s1.SHRTRCE_SUBJ_CODE = 'MATH'
AND s1.SHRTRCE_CRSE_NUMB = '2550'
AND s1.SHRTRCE_CREDIT_HOURS > 0
) ON M2550PIDM = pop.PIDM

--checks for MATH2551 pre-enrollment credit
LEFT OUTER JOIN
(
SELECT DISTINCT 
s1.SHRTRCE_PIDM AS M2551PIDM,
s1.SHRTRCE_TERM_CODE_EFF AS M2551TERM
FROM EDWSTG.STG_BAN_SHRTRCE_TBL s1
INNER JOIN EDWSTG.STG_BAN_SHRTRCR_TBL s2
ON s1.SHRTRCE_PIDM = s2.SHRTRCR_PIDM AND s1.SHRTRCE_TRCR_SEQ_NO = s2.SHRTRCR_SEQ_NO
INNER JOIN EDWSTG.STG_BAN_SHRTRIT_TBL s3
ON s3.SHRTRIT_PIDM = s2.SHRTRCR_PIDM AND s2.SHRTRCR_TRIT_SEQ_NO = s3.SHRTRIT_SEQ_NO
WHERE s1.SHRTRCE_SUBJ_CODE = 'MATH'
AND s1.SHRTRCE_CRSE_NUMB = '2551'
AND s1.SHRTRCE_CREDIT_HOURS > 0
) ON M2551PIDM = pop.PIDM

--checks for CS1301 pre-enrollment credit
LEFT OUTER JOIN
(
SELECT DISTINCT 
s1.SHRTRCE_PIDM AS CS1301PIDM,
s1.SHRTRCE_TERM_CODE_EFF AS CS1301TERM
FROM EDWSTG.STG_BAN_SHRTRCE_TBL s1
INNER JOIN EDWSTG.STG_BAN_SHRTRCR_TBL s2
ON s1.SHRTRCE_PIDM = s2.SHRTRCR_PIDM AND s1.SHRTRCE_TRCR_SEQ_NO = s2.SHRTRCR_SEQ_NO
INNER JOIN EDWSTG.STG_BAN_SHRTRIT_TBL s3
ON s3.SHRTRIT_PIDM = s2.SHRTRCR_PIDM AND s2.SHRTRCR_TRIT_SEQ_NO = s3.SHRTRIT_SEQ_NO
WHERE s1.SHRTRCE_SUBJ_CODE = 'CS'
AND s1.SHRTRCE_CRSE_NUMB = '1301'
AND s1.SHRTRCE_CREDIT_HOURS > 0
) ON CS1301PIDM = pop.PIDM

--checks for CS1331 pre-enrollment credit
LEFT OUTER JOIN
(
SELECT DISTINCT 
s1.SHRTRCE_PIDM AS CS1331PIDM,
s1.SHRTRCE_TERM_CODE_EFF AS CS1331TERM
FROM EDWSTG.STG_BAN_SHRTRCE_TBL s1
INNER JOIN EDWSTG.STG_BAN_SHRTRCR_TBL s2
ON s1.SHRTRCE_PIDM = s2.SHRTRCR_PIDM AND s1.SHRTRCE_TRCR_SEQ_NO = s2.SHRTRCR_SEQ_NO
INNER JOIN EDWSTG.STG_BAN_SHRTRIT_TBL s3
ON s3.SHRTRIT_PIDM = s2.SHRTRCR_PIDM AND s2.SHRTRCR_TRIT_SEQ_NO = s3.SHRTRIT_SEQ_NO
WHERE s1.SHRTRCE_SUBJ_CODE = 'CS'
AND s1.SHRTRCE_CRSE_NUMB = '1331'
AND s1.SHRTRCE_CREDIT_HOURS > 0
) ON CS1331PIDM = pop.PIDM

--checks for CS1371 pre-enrollment credit
LEFT OUTER JOIN
(
SELECT DISTINCT 
s1.SHRTRCE_PIDM AS CS1371PIDM,
s1.SHRTRCE_TERM_CODE_EFF AS CS1371TERM
FROM EDWSTG.STG_BAN_SHRTRCE_TBL s1
INNER JOIN EDWSTG.STG_BAN_SHRTRCR_TBL s2
ON s1.SHRTRCE_PIDM = s2.SHRTRCR_PIDM AND s1.SHRTRCE_TRCR_SEQ_NO = s2.SHRTRCR_SEQ_NO
INNER JOIN EDWSTG.STG_BAN_SHRTRIT_TBL s3
ON s3.SHRTRIT_PIDM = s2.SHRTRCR_PIDM AND s2.SHRTRCR_TRIT_SEQ_NO = s3.SHRTRIT_SEQ_NO
WHERE s1.SHRTRCE_SUBJ_CODE = 'CS'
AND s1.SHRTRCE_CRSE_NUMB = '1371'
AND s1.SHRTRCE_CREDIT_HOURS > 0
) ON CS1371PIDM = pop.PIDM

--checks for PHYS2211 pre-enrollment credit
LEFT OUTER JOIN
(
SELECT DISTINCT 
s1.SHRTRCE_PIDM AS PHYS2211PIDM,
s1.SHRTRCE_TERM_CODE_EFF AS PHYS2211TERM
FROM EDWSTG.STG_BAN_SHRTRCE_TBL s1
INNER JOIN EDWSTG.STG_BAN_SHRTRCR_TBL s2
ON s1.SHRTRCE_PIDM = s2.SHRTRCR_PIDM AND s1.SHRTRCE_TRCR_SEQ_NO = s2.SHRTRCR_SEQ_NO
INNER JOIN EDWSTG.STG_BAN_SHRTRIT_TBL s3
ON s3.SHRTRIT_PIDM = s2.SHRTRCR_PIDM AND s2.SHRTRCR_TRIT_SEQ_NO = s3.SHRTRIT_SEQ_NO
WHERE s1.SHRTRCE_SUBJ_CODE = 'PHYS'
AND s1.SHRTRCE_CRSE_NUMB = '2211'
AND s1.SHRTRCE_CREDIT_HOURS > 0
) ON PHYS2211PIDM = pop.PIDM

--checks for PHYS2212 pre-enrollment credit
LEFT OUTER JOIN
(
SELECT DISTINCT 
s1.SHRTRCE_PIDM AS PHYS2212PIDM,
s1.SHRTRCE_TERM_CODE_EFF AS PHYS2212TERM
FROM EDWSTG.STG_BAN_SHRTRCE_TBL s1
INNER JOIN EDWSTG.STG_BAN_SHRTRCR_TBL s2
ON s1.SHRTRCE_PIDM = s2.SHRTRCR_PIDM AND s1.SHRTRCE_TRCR_SEQ_NO = s2.SHRTRCR_SEQ_NO
INNER JOIN EDWSTG.STG_BAN_SHRTRIT_TBL s3
ON s3.SHRTRIT_PIDM = s2.SHRTRCR_PIDM AND s2.SHRTRCR_TRIT_SEQ_NO = s3.SHRTRIT_SEQ_NO
WHERE s1.SHRTRCE_SUBJ_CODE = 'PHYS'
AND s1.SHRTRCE_CRSE_NUMB = '2212'
AND s1.SHRTRCE_CREDIT_HOURS > 0
) ON PHYS2212PIDM = pop.PIDM

--check graduation status
LEFT OUTER JOIN
(
SELECT
o2.PIDM AS GRADPIDM,
o2.GRADUATION_TERM_CODE AS GRADTERM,
o2.MAJOR AS GRAD_Major,
o2.DEGR_CODE AS GRAD_Awarded
FROM (
       SELECT o.*, row_number() OVER (PARTITION BY o.PIDM ORDER BY o.GRADUATION_TERM_CODE ASC) AS rn
       FROM EDW.AGG_STU_OUTCOME_MVW o WHERE o.LVL_CODE = 'US' AND o.DEGR_CODE LIKE 'B%' AND o.GRADUATION_TERM_CODE > 201408
     )o2
WHERE o2.rn = 1
) ON GRADPIDM = pop.PIDM

WHERE FIRST_TERM IS NOT NULL
AND FIRST_TERM >= '201405'
AND HRS_EARNED BETWEEN 1 AND 10000