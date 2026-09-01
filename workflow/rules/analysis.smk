"""Ordered analysis stages executed through the project pipeline."""


DETERMINISTIC_ENV = (
    f"PYTHONHASHSEED={RANDOM_SEED} "
    "LANG=C.UTF-8 "
    "LC_ALL=C.UTF-8 "
    "TZ=UTC "
    "SOURCE_DATE_EPOCH=1788220800 "
    "OMP_NUM_THREADS=1 "
    "OPENBLAS_NUM_THREADS=1 "
    "MKL_NUM_THREADS=1 "
    "NUMEXPR_NUM_THREADS=1 "
    "MPLBACKEND=Agg"
)
rule validate:
    input:
        config=CONFIG_PATH,
        manifest=MANIFEST_PATH,
        local_inputs=LOCAL_INPUTS,
        schemas=SCHEMA_INPUTS,
        analysis_sources=ANALYSIS_SOURCES,
        workflow_sources=WORKFLOW_SOURCES,
        environment=ENVIRONMENT_INPUTS,
    output:
        marker=ensure(str(Path(STATUS_DIR) / "validate.done"), non_empty=True),
        report=ensure(VALIDATION_OUTPUT, non_empty=True),
    log:
        str(Path(OUTPUT_DIR) / "logs" / "validate.log")
    params:
        env=DETERMINISTIC_ENV,
        output_dir=lambda _wildcards, output: str(Path(output[0]).parents[1]),
        python=PYTHON_EXECUTABLE,
        random_seed=RANDOM_SEED,
    threads: 1
    shell:
        """
        {params.env} {params.python:q} -m chemokine_cart.pipeline \
            --stage validate \
            --config {input.config:q} \
            --manifest {input.manifest:q} \
            --random-seed {params.random_seed} \
            --output-dir {params.output_dir:q} > {log:q} 2>&1
        test -s {output.marker:q}
        """


rule single_cell:
    input:
        previous=rules.validate.output.marker,
        config=CONFIG_PATH,
        manifest=MANIFEST_PATH,
    output:
        marker=ensure(str(Path(STATUS_DIR) / "single_cell.done"), non_empty=True),
        tables=ensure(SINGLE_CELL_OUTPUTS, non_empty=True),
    log:
        str(Path(OUTPUT_DIR) / "logs" / "single_cell.log")
    params:
        env=DETERMINISTIC_ENV,
        output_dir=lambda _wildcards, output: str(Path(output[0]).parents[1]),
        python=PYTHON_EXECUTABLE,
        random_seed=RANDOM_SEED,
    threads: 1
    shell:
        """
        {params.env} {params.python:q} -m chemokine_cart.pipeline \
            --stage single-cell \
            --config {input.config:q} \
            --manifest {input.manifest:q} \
            --random-seed {params.random_seed} \
            --output-dir {params.output_dir:q} > {log:q} 2>&1
        test -s {output.marker:q}
        """


rule spatial:
    input:
        previous=rules.single_cell.output.marker,
        config=CONFIG_PATH,
        manifest=MANIFEST_PATH,
    output:
        marker=ensure(str(Path(STATUS_DIR) / "spatial.done"), non_empty=True),
        tables=ensure(SPATIAL_OUTPUTS, non_empty=True),
    log:
        str(Path(OUTPUT_DIR) / "logs" / "spatial.log")
    params:
        env=DETERMINISTIC_ENV,
        output_dir=lambda _wildcards, output: str(Path(output[0]).parents[1]),
        python=PYTHON_EXECUTABLE,
        random_seed=RANDOM_SEED,
    threads: 1
    shell:
        """
        {params.env} {params.python:q} -m chemokine_cart.pipeline \
            --stage spatial \
            --config {input.config:q} \
            --manifest {input.manifest:q} \
            --random-seed {params.random_seed} \
            --output-dir {params.output_dir:q} > {log:q} 2>&1
        test -s {output.marker:q}
        """


rule functional:
    input:
        previous=rules.spatial.output.marker,
        config=CONFIG_PATH,
        manifest=MANIFEST_PATH,
    output:
        marker=ensure(str(Path(STATUS_DIR) / "functional.done"), non_empty=True),
        tables=ensure(FUNCTIONAL_OUTPUTS, non_empty=True),
    log:
        str(Path(OUTPUT_DIR) / "logs" / "functional.log")
    params:
        env=DETERMINISTIC_ENV,
        output_dir=lambda _wildcards, output: str(Path(output[0]).parents[1]),
        python=PYTHON_EXECUTABLE,
        random_seed=RANDOM_SEED,
    threads: 1
    shell:
        """
        {params.env} {params.python:q} -m chemokine_cart.pipeline \
            --stage functional \
            --config {input.config:q} \
            --manifest {input.manifest:q} \
            --random-seed {params.random_seed} \
            --output-dir {params.output_dir:q} > {log:q} 2>&1
        test -s {output.marker:q}
        """


rule figures:
    input:
        previous=rules.functional.output.marker,
        config=CONFIG_PATH,
        manifest=MANIFEST_PATH,
    output:
        marker=ensure(str(Path(STATUS_DIR) / "figures.done"), non_empty=True),
        figures=ensure(FIGURE_OUTPUTS, non_empty=True),
    log:
        str(Path(OUTPUT_DIR) / "logs" / "figures.log")
    params:
        env=DETERMINISTIC_ENV,
        output_dir=lambda _wildcards, output: str(Path(output[0]).parents[1]),
        python=PYTHON_EXECUTABLE,
        random_seed=RANDOM_SEED,
    threads: 1
    shell:
        """
        {params.env} {params.python:q} -m chemokine_cart.pipeline \
            --stage figures \
            --config {input.config:q} \
            --manifest {input.manifest:q} \
            --random-seed {params.random_seed} \
            --output-dir {params.output_dir:q} > {log:q} 2>&1
        test -s {output.marker:q}
        """


rule report:
    input:
        previous=rules.figures.output.marker,
        config=CONFIG_PATH,
        manifest=MANIFEST_PATH,
    output:
        marker=ensure(str(Path(STATUS_DIR) / "report.done"), non_empty=True),
        reports=ensure(REPORT_OUTPUTS, non_empty=True),
    log:
        str(Path(OUTPUT_DIR) / "logs" / "report.log")
    params:
        env=DETERMINISTIC_ENV,
        output_dir=lambda _wildcards, output: str(Path(output[0]).parents[1]),
        python=PYTHON_EXECUTABLE,
        random_seed=RANDOM_SEED,
    threads: 1
    shell:
        """
        {params.env} {params.python:q} -m chemokine_cart.pipeline \
            --stage report \
            --config {input.config:q} \
            --manifest {input.manifest:q} \
            --random-seed {params.random_seed} \
            --output-dir {params.output_dir:q} > {log:q} 2>&1
        test -s {output.marker:q}
        """
