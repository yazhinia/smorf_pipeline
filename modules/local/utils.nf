// Nextflow's task cache does not see edits to bin/ scripts. Each Python step takes the
// scripts it runs as an extra `code` input (staged, never read from there), so editing a
// helper re-runs exactly the tasks that use it on -resume.
def binFiles(List names) {
    return names.collect { n -> file("${projectDir}/bin/${n}", checkIfExists: true) }
}
