#!/usr/bin/env python3
# \HEADER\-------------------------------------------------------------------------
#
#  CONTENTS      : snakemake mxq submission script
#
#  DESCRIPTION   : none
#
#  RESTRICTIONS  : none
#
#  REQUIRES      : none
#
# ---------------------------------------------------------------------------------
# Copyright (c) 2018-2019, Pay Giesselmann, Max Planck Institute for Molecular Genetics
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# Written by Pay Giesselmann, changed by Anton Kriese
# ---------------------------------------------------------------------------------
import sys, subprocess
from snakemake.utils import read_job_properties

if __name__ == '__main__':
    jobscript = sys.argv[-1]
    job_properties = read_job_properties(jobscript)

    # default resources
    threads = '10'
    runtime = '600'
    memory = '10000'
    tmpdir = '10' # default of mxq are 10G
    # parse resources
    if job_properties["type"] == "group":
        group_name = job_properties["groupid"]
    else:
        group_name = job_properties["rule"]
    if "threads" in job_properties:
        threads = str(job_properties["threads"])
    resources = job_properties["resources"]
    if "time_min" in resources:
        runtime = str(resources["time_min"])
    if "mem_mb" in resources:
        memory = str(resources["mem_mb"])
    if "tmpdir_gb" in resources:
        tmpdir = str(resources["tmpdir_gb"])
        # print(tmpdir, file=sys.stderr)

    # cmd and submission
    blacklist = ' '.join(['dontpanic', 'fordprefect'])
    cmd = f'mxqsub --blacklist "{blacklist}" --threads={threads} --memory={memory} --runtime={runtime} --tmpdir {tmpdir}G --group-name={group_name} --command-alias={group_name} {jobscript}'
    try:
        ret = subprocess.run(cmd, check=True, shell=True, stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        raise e

    # parse job ID from output
    out = ret.stdout.decode()
    sub = {key:value for key, value in [line.split('=') for line in out.strip().split()]}

    # print job id for snakemake
    print(sub['mxq_job_id'])
