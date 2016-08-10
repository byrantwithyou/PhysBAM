;#####################################################################
; Copyright 2004-2006.
;#####################################################################
; physbam.el
;#####################################################################

;#####################################################################
; PhysBAM Style
;#####################################################################

(defconst physbam-c-style
  '((c-basic-offset . 4)
        (c-offsets-alist . ((string . -1000)
                        (c . c-lineup-C-comments)
                        (defun-open . 0)
                        (defun-close . 0)
                        (defun-block-intro . +)
                        (class-open . 0)
                        (class-close . 0)
                        (inline-open . 0)
                        (inline-close . 0)
                        (func-decl-cont . +)
                        (knr-argdecl-intro . 5)
                        (knr-argdecl . 0)
                        (topmost-intro . 0)
                        (topmost-intro-cont . 0)
                        (member-init-intro . +)
                        (member-init-cont . -1)
                        (inher-intro . +)
                        (inher-cont . c-lineup-multi-inher)
                        (block-open . 0)
                        (block-close . 0)
                        (brace-list-open . 0)
                        (brace-list-close . 0)
                        (brace-list-intro . +)
                        (brace-list-entry . 0)
                        (statement . 0)
                        (statement-cont . +)
                        (statement-block-intro . +)
                        (statement-case-intro . +)
                        (statement-case-open . +)
                        (substatement . +)
                        (substatement-open . 0)
                        (case-label . +)
                        (access-label . -)
                        (label . *)
                        (do-while-closure . 0)
                        (else-clause . 0)
                        (comment-intro . c-lineup-comment)
                        (arglist-intro . +)
                        (arglist-cont . 0)
                        (arglist-cont-nonempty . +)
                        (arglist-close . 0)
                        (stream-op . c-lineup-streamop)
                        (inclass . +)
                        (cpp-macro . -1000)
                        (friend . 0)
                        (extern-lang-open . 0)
                        (extern-lang-close . 0)
                        (inextern-lang . +)
			(namespace-open . 0)
			(namespace-close . 0)
			(innamespace . 0)
                        (template-args-cont . +)))    
    (c-hanging-braces-alist . ((class-open before after)
                               (class-close before)
                               (defun-open before after)
                               (defun-close before after)
                               (inline-open before)
                               (inline-close after)
                               (brace-list-open)
                               (brace-list-close)
                               (brace-list-intro)
                               (brace-list-entry)
                               (block-open after)
                               (block-close)
                               (substatement-open after)
                               (substatement-case-open)
			       (namespace-open after)
			       (namespace-close before after)
                               (extern-lang-open after)
                               (extern-lang-close after)))
    (c-hanging-colons-alist . ((case-label)
                               (label)
                               (access-label after)
                               (member-init-intro before)
                               (inher-intro)))
    (c-cleanup-list         . ((defun-close-semi)
                               (scope-operator)))
    (c-hanging-semi&comma-criteria . (c-semi&comma-no-newlines-before-nonblanks))
    (c-echo-syntactic-information-p . t))
  "PhysBAM C++ Style")

(c-add-style "physbam" physbam-c-style)

;#####################################################################
; PhysBAM Helper Routines
;#####################################################################

(defun physbam-reduce (f x)
  (if (eq (cdr x) nil)
      (car x)
    (funcall f (car x) (physbam-reduce f (cdr x))))) 

(defun physbam-filter (f list)
  (let (filtered)
    (dolist (x list filtered)
      (if (funcall f x) (setq filtered (cons x filtered)) nil))))

;#####################################################################
; PhysBAM Navigation Commands
;#####################################################################

(defun physbam-header-flip ()
  "Find the .h file for this .C file (or vice versa)."
  (interactive)
  (let ((dotc (string-match "[.]\\(cpp\\|c\\)$" (buffer-file-name)))
        (doth (string-match "[.]h$" (buffer-file-name))))
    (if dotc
        (find-file (concat (substring (buffer-file-name) 0 dotc) ".h"))
      (if doth
          (find-file (concat (substring (buffer-file-name) 0 doth) ".cpp"))
        (message "Not a cpp or h file!!")))))

(defun physbam-open-parent ()
  "Open header of parent class"
  (interactive)
  (let ((doth (string-match "[.]h$" (buffer-file-name))))
    (if doth
        (let ((save_point (point)))
          (goto-char 0)
          (if (re-search-forward ":\\(?:public\\|protected\\|private\\) \\(\\(?:\\w\\|_\\)*\\)" nil t)
              (let ((filename (concat (match-string 1) ".h")))
                (goto-char save_point)
                (if (file-exists-p filename)
                    (find-file filename)
                  (let ((physbam_filename (find-physbam-file filename)))
                    (if (file-exists-p physbam_filename)
                        (find-file physbam_filename)
                      (message "Can't find header of parent")))))
            (goto-char save_point)
            (message "No parent class")))
      (message "Not a header file"))))

(defun physbam-find-class ()
  "Open header of parent class"
  (interactive)
  (let ((save_point (point)))
    (search-forward-regexp "[^A-Za-z0-9_]")
    (backward-char)
    (let ((a (point)))
      (search-backward-regexp "[^A-Za-z0-9_]")
      (forward-char)
      (let ((b (point)))
        (let ((class-name (buffer-substring a b)))
          (let ((filename (concat class-name ".h")))
            (goto-char save_point)
            (if (file-exists-p filename)
                (find-file filename)
              (let ((physbam-filename (find-physbam-file filename)))
                (if (file-exists-p physbam-filename)
                    (find-file physbam-filename)
                  (message (concat "Can't find " filename)))))))))))

(defun find-physbam-file (filename)
  (call-process "find" nil "*physbam-filename*" nil  (getenv "PUBLIC") "-name" filename)
  (set-buffer "*physbam-filename*")
  (goto-char 0)
  (let ((found_filename (buffer-substring (point) (line-end-position))))
    (kill-buffer "*physbam-filename*")
    (if (file-exists-p found_filename) found_filename nil)))

(defun physbam-grep-cpp-and-headers (directory_prefix querystr)
  (let ((path (concat (getenv "PHYSBAM") "/" directory_prefix "/")))
  (grep (concat "(find " path " -name '*.h' -o -name '*.hpp' -o -name '*.cpp' | xargs grep -n -e \"" querystr "\") # "))))

(defun physbam-grep-public (querystr)
  "grep in Public_Library"
  (interactive "sGrep $PHYSBAM/Public_Library for: ")
  (physbam-grep-cpp-and-headers "Public_Library" querystr))

(defun physbam-grep-projects (querystr)
  "grep in Projects"
  (interactive "sGrep $PHYSBAM/Projects for:")
  (physbam-grep-cpp-and-headers "*Projects" querystr))

(defun physbam-grep-tools (querystr)
  "grep in Tools"
  (interactive "sGrep $PHYSBAM/Tools for:")
  (physbam-grep-cpp-and-headers "Tools" querystr))

(defun physbam-fix-includes ()
  "Fix includes in current file"
  (interactive)
  (save-buffer)
  (shell-command (format "%s/Scripts/misc/fix_headers_file.sh %s" (getenv "PHYSBAM") (buffer-file-name)))
  (revert-buffer-no-prompt))

(defun grep-cpp-and-headers (querystr)
  "grep headers and source files recursively in current directory"
  (interactive "sGrep for:")
  (grep (concat "(find . -name '*.h' -o -name '*.hpp' -o -name '*.cpp' | xargs grep -n -e \"" querystr "\") # ")))

;#####################################################################
; PhysBAM Formatting Helpers
;#####################################################################

(defun physbam-fix-function-comment ()
  "Fix the physbam style function comments"
  (interactive)
  (setq function-name (physbam-get-function-name))
  (setq function-type (physbam-get-function-type function-name))
  (let ((line-count 0))
    (beginning-of-line)
    (if (looking-at "[a-zA-Z0-9_]") 
        (progn (forward-line -1)
               (setq line-count (+ line-count 1))))
    (while (looking-at "//") 
      (if (looking-at "//#")
          (while (looking-at "//")
            (delete-region (point) (line-end-position))
            (delete-char 1)
            (forward-line -1))
        (progn
          (setq line-count (+ line-count 1))
          (forward-line -1))))
    (forward-line +1)
    (physbam-insert-function-comment function-name function-type)
    (forward-line line-count)))

(defun physbam-insert-function-comment (name type)
  "Inserts a physbam style function comment at point"
  (interactive)
  (insert "//#####################################################################\n")
  (insert (concat "// " type))
  (if (string= type "Function")
      (insert (concat " " name)))
  (insert "\n")
  (insert "//#####################################################################\n"))

(defun physbam-get-function-name ()
  "Gets the name of the function on the current line"
  (beginning-of-line)
  (let ((curr-line (buffer-substring (point) (line-end-position))))
    (string-match "^\\(~?[a-zA-Z0-9_]+\\)\(.+$" curr-line)
    (match-string 1 curr-line)))

(defun physbam-get-function-type (name)
  "Determines the type of function on the current line"
  (forward-line -1)
  (beginning-of-line)
  (let ((curr-line (buffer-substring (point) (line-end-position))))
    (cond ((string-match "~" name) "Destructor")
          ((string-match name curr-line) "Constructor")
          (t "Function"))))

(defun physbam-get-lastname (fullname)
  (car (last (split-string fullname "[ ]+"))))

(defun physbam-current-year ()
  "Get year"
  (string-to-number (substring (current-time-string) -4)))

(defun physbam-fix-copyright-user (fullname)
  "Fix the list of names in copyright i.e. insert mine"
  (interactive)
  (goto-line 2)
  (beginning-of-line)
  (let ((first-line (buffer-substring (point) (line-end-position))))
    (let ((items (split-string first-line "[ \t]*,[ \t]*")))
      (let ((names '()) (years '()) (year-regexp "\\([0-9]\\{4\\}\\)"))
        (dolist (x items)
          (if (string-match (format "^\\(// Copyright[ ]\\)?+%s$" year-regexp) x)
              (setq years (cons (match-string 2 x) years))
            (if (string-match (format "^\\(// Copyright[ ]\\)?+%s-%s$" year-regexp year-regexp) x)
                (setq years (cons (match-string 2 x) (cons (match-string 3 x) years)))
              (setq names (cons (replace-regexp-in-string  "\\." "" x) names)))))
        (unless (member fullname names) (setq names (cons fullname names)))
        (let ((sorted-names (sort names '(lambda (x y) (string< (physbam-get-lastname x) (physbam-get-lastname y))))))
          (delete-region (point) (line-end-position))
          (insert (format "// Copyright %s, %s." 
                          (let ((min-year (physbam-reduce (lambda (x y) (if (< x y) x y)) (mapcar 'string-to-number years)))
                                (max-year (physbam-current-year)))
                            (if (= min-year max-year) (format "%d" min-year) (format "%d-%d" min-year max-year)))
                          (physbam-reduce (lambda (x y) (concat x ", " y)) sorted-names))))))))

(defun physbam-fix-copyright ()
  (interactive)
  (physbam-fix-copyright-user user-full-name))

(defun physbam-fix-copyright-user-list ()
  (interactive)
  (dolist (fullname user-list)
    (physbam-fix-copyright-user fullname)))

(defun physbam-insert-copyright ()
  (interactive)
  (goto-char 0) (goto-line 7)
  (insert "//#####################################################################\n")
  (insert (concat "// Copyright " (number-to-string (physbam-current-year))  ", " user-full-name ".\n"))
  (insert "// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.\n")
  (insert "//#####################################################################\n"))

(defun physbam-insert-header ()
  (interactive)
  (physbam-insert-copyright)
  (let ((classname (substring (buffer-name) 0 -2)))
    (insert (format "#ifndef __%s__\n" classname))
    (insert (format "#define __%s__\n\n" classname))
    (insert "namespace PhysBAM{\n\n")
    (insert "template<class T>\n")
    (insert (format "class %s\n" classname))
    (insert "{\n")
    (insert "public:\n\n")
    (insert "//#####################################################################\n")
    (insert "};\n")
    (insert "}\n")
    (insert "#endif\n")))

(defun physbam-add-include ()
  "Add include from current location"
  (interactive)
  (search-forward-regexp "[^A-Za-z0-9_]")
  (backward-char)
  (setq a (point))
  (search-backward-regexp "[^A-Za-z0-9_]")
  (forward-char)
  (setq b (point))
  (setq class-name (buffer-substring a b))
  (goto-line 0)
  (beginning-of-line)
  (search-forward "#include")
  (beginning-of-line)
  (insert (concat "#include <temp/" class-name ".h>\n"))
  (physbam-fix-includes))

(defun physbam-make-cpp-from-header ()
  (interactive)
  (goto-char 0)
  (goto-line 7)
  (copy-region-as-kill (point-min) (point))
  (physbam-header-flip)
  (yank)
  (let ((classname (substring (buffer-name) 0 -4)))
    (insert (format "#include \"%s.h\"\n" classname))
    (insert "using namespace PhysBAM;\n")
    (insert (format "template class %s<float>;\n" classname))
    (insert (format "template class %s<double>;\n" classname))))

;#####################################################################
; PhysBAM Build / Run Commands
;#####################################################################
; returns the last eleement of the path i.e. the project

(defun physbam-set-project-type (type)
  (setq physbam-project-type type)
  (physbam-setup-compile-command))

(defun physbam-set-project-directory (project-directory)
  (setq physbam-project-directory project-directory)
  (setq physbam-project-directory (concat physbam-project-directory (if (string= (substring physbam-project-directory -1 nil) "/") "" "/")))
  (physbam-read-project-settings)
  (physbam-setup-compile-command))

(defun physbam-choose-project-directory ()
  (interactive)
  (physbam-set-project-directory (read-file-name "New Project Directory: " physbam-project-directory () nil)))
  

(defun physbam-choose-output-directory ()
  (interactive)
  (setq physbam-output-directory (read-file-name "New output directory: " physbam-output-directory () nil))
  (physbam-setup-compile-command))

(defun physbam-set-compile-count (count)
  (interactive)
  (setq physbam-compile-count count)
  (physbam-setup-compile-command))

(defun physbam-compile ()
  (interactive)
  (let ((old-default-directory default-directory))
    (setq default-directory physbam-project-directory)
    (call-interactively 'compile)
    (setq default-directory old-default-directory)))

(setq compile-current-file-command "g++ -c -g -Wno-unused-local-typedefs -Wall -Werror -Winit-self -Woverloaded-virtual -Wstrict-aliasing=2 -std=gnu++14 -Wno-unknown-pragmas -Wno-strict-overflow -Wno-sign-compare -I$PHYSBAM/Public_Library -o /dev/null ")
(defun physbam-compile-current-file ()
  (interactive)
  (let ((old-compile-command compile-command))
    (setq compile-command (concat compile-current-file-command (buffer-file-name)))
    (save-buffer)
    (call-interactively 'compile)
    (setq compile-command old-compile-command)))

(defun physbam-setup-compile-command ()
  (setq compile-command (format "if [ -e Makefile ] ; then make -k -j %d ; else nice scons --warn=no-duplicate-environment --warn=no-deprecated -Q --implicit-cache -k -u TYPE=%s -j %d ; fi" physbam-compile-count physbam-project-type physbam-compile-count))
  (message (format "New compile command is: %s" compile-command))
 (define-key global-map [menu-bar physbam project-dir] '(menu-item (concat  "Project: " physbam-project-directory) physbam-choose-project-directory)))

;#####################################################################
; Menu Bar
;#####################################################################
; Generic helpful scripts
(define-key global-map [menu-bar physbam] (cons "PhysBAM" (make-sparse-keymap "PhysBAM")))
(define-key global-map [menu-bar physbam edit-physbam-settings] '("Edit PhysBAM Emacs" . (lambda () (interactive) (find-file (concat (getenv "PHYSBAM") "/Scripts/emacs/physbam.el")))))
(define-key global-map [menu-bar physbam edit-settings] '("Edit .emacs" . (lambda () (interactive) (find-file (concat (getenv "PHYSBAM") "~/.emacs")))))
(define-key global-map [menu-bar physbam sep0] '(menu-item "--single-line"))
(define-key global-map [menu-bar physbam fix-copyright] '("Fix Copyright" . physbam-fix-copyright))
(define-key global-map [menu-bar physbam fix-function-comment] '("Fix Function Comment" . physbam-fix-function-comment))
(define-key global-map [menu-bar physbam sep1] '(menu-item "--single-line"))
; other compile stuff
(define-key global-map [menu-bar physbam compile] '("Compile Project" . physbam-compile))
(define-key global-map [menu-bar physbam compile-file] '("Compile Current Buffer" . physbam-compile-current-file))
(define-key global-map [menu-bar physbam sep1b] '(menu-item "--single-line"))
; Debug or release
(define-key global-map [menu-bar physbam release] '(menu-item "Release" (lambda () (interactive) (physbam-set-project-type "release")) :button (:toggle . (string= physbam-project-type "release"))))
(define-key global-map [menu-bar physbam optdebug] '(menu-item "Optimized Debug" (lambda () (interactive) (physbam-set-project-type "optdebug")) :button (:toggle . (string= physbam-project-type "optdebug"))))
(define-key global-map [menu-bar physbam debug] '(menu-item "Debug" (lambda () (interactive) (physbam-set-project-type "debug")) :button (:toggle . (string= physbam-project-type "debug"))))
(define-key global-map [menu-bar physbam profile] '(menu-item "Profile" (lambda () (interactive) (physbam-set-project-type "profile")) :button (:toggle . (string= physbam-project-type "profile"))))
(define-key global-map [menu-bar physbam sep2] '(menu-item "--single-line"))
(define-key global-map [menu-bar physbam autosetproject] '(menu-item "Set Project From Current" (lambda () (interactive) (physbam-find-project))))

;#####################################################################
; Development related variables
;#####################################################################

(defun physbam-find-project-helper (filename orig-filename depth)
  (cond
   ((file-exists-p (concat filename "SConscript")) (expand-file-name filename))
   ((> depth 10) orig-filename)
   (t (physbam-find-project-helper (concat filename "../") orig-filename (+ depth 1)))))

(defun physbam-find-project ()
  (interactive)
  (setq physbam-project-directory (physbam-find-project-helper default-directory default-directory 0)))

(physbam-find-project)

(setq physbam-project-type "release")
(setq physbam-compile-count (let ((count (getenv "PHYSBAM_COMPILE_COUNT"))) (if count (string-to-number count) 4)))
(physbam-setup-compile-command)
(setq truncate-partial-width-windows nil)
(setq compilation-scroll-output t) ; scroll to end by default
(setq compilation-read-command nil)
