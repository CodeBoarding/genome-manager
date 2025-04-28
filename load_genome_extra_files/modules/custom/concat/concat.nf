/*
Concatenate two file and merge metadata
*/

process CONCAT {
    label 'process_1cpu_512mb_30min'

    tag "${meta1.id}:${meta2.id}"

    input:
    tuple val(meta1), path(file1)
    tuple val(meta2), path(file2)

    output:
    tuple val(meta), path("${file1.toString().take(file1.toString().lastIndexOf('.'))}.${file2.toString()}"), emit: merged

    script:
    meta = meta1.clone()
    meta.id = meta1.id + ':' + meta2.id
    meta.id1 = meta1.id
    meta.id2 = meta2.id
    // def file1_name = file1.toString()
    // def file2_name = file2.toString()
    // outname = file1_name.take(file1_name.lastIndexOf('.')) + '.' + file2_name
    """
    cat $file1 $file2 > ${file1.toString().take(file1.toString().lastIndexOf('.'))}.${file2.toString()}
    """
}
