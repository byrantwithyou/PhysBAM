VIM PhysBAM Extensions

* Install
  Require `git` and `grep` installed.
  Source `physbam.vim` in your `.vimrc`.
  Bind functions to your prefered keys.

  Example:
  
    nnoremap <F2> :call OpenFilePair()<CR>
    nnoremap <F3> :call OpenClassHeader()<CR>
    nnoremap <F4> :call CopyClassHeader()<CR>

* Functions list
+ `OpenFilePair`
  Switch between header/cpp pair.
+ `OpenClassHeader`
  Open the header for the class name under the cursor.
+ `CopyClassHeader`
  Generate a cpp include line (i.e. "#include ...") and store it in `h` register.

