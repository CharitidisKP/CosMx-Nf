// Samplesheet rows with include = FALSE are skipped everywhere. A missing column counts as TRUE.
def isIncluded(row) {
    def v = (row.include ?: 'TRUE').toString().trim().toUpperCase()
    if (!(v in ['TRUE', 'FALSE']))
        error "samplesheet: include must be TRUE or FALSE, got '${row.include}' for ${row.sample_id}"
    return v == 'TRUE'
}
