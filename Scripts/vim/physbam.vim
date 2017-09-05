function! FilePair()
    if expand('%:e') ==# 'cpp'
        return expand('%:p:h') . '/' . expand('%:r') . '.h'
    elseif expand('%:e') ==# 'h'
        return expand('%:p:h') . '/' . expand('%:r') . '.cpp'
    else
        return ''
    endif
endfunction

function! TryOpenFile(filename)
    if filereadable(a:filename)
        echom 'open ' . a:filename
        execute 'edit' a:filename
    endif
endfunction

function! OpenFilePair()
    let filename = FilePair()
    if filename !=# ''
        call TryOpenFile(filename)
    endif
endfunction

function! ClassName()
    return expand('<cword>')
endfunction

function! ClassHeader()
    let root = systemlist('git rev-parse --show-toplevel')[0]
    let name = '__' . ClassName() . '__'
    let cmd = 'grep -nrl ' . name . ' ' . root . ' --include \*.h'
    let matches = systemlist(cmd)
    if len(matches) == 1
        return matches[0]
    else
        return ''
    endif
endfunction

function! OpenClassHeader()
    let filename = ClassHeader()
    if filename !=# ''
        call TryOpenFile(filename)
    endif
endfunction

function! CopyClassHeader()
    let filename = ClassHeader()
    if filename !=# ''
        let init = systemlist('git rev-parse --show-toplevel')[0] . '/Public_Library/'
        let header = substitute(filename, init, '', 'g')
        echom header
        let @h = '#include <' . header . ">\n"
    endif
endfunction

