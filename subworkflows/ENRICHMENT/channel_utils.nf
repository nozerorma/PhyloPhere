#!/usr/bin/env nextflow
// channel_utils.nf — helpers shared by the FCS sub-workflows.

// A channel that several processes read must be a value channel. A value channel is returned as it is (first() on
// it is redundant and Nextflow warns about it); a queue channel is reduced to its first item.
def asValue(ch) {
    return (ch instanceof groovyx.gpars.dataflow.DataflowVariable) ? ch : ch.first()
}
