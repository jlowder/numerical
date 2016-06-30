(in-package :cl-user)

(defpackage :numerical
  (:use :common-lisp)
  (:export :if_
           :when_
           :cond_
           :_
           :@))

(in-package :numerical)

;; define some constants to be used as defaults throughout the library. 
(defconstant +limit+ 50) ; limit on the number of attempts
(defconstant +tolerance+ 5d-8) ; approximate to 8 decimal places of precision
(defconstant +epsilon+ 5d-20) ; amount to nudge to make non-zero

;; define some anaphoric macros that can optionally be used with
;; functions throughout this library, which return 3 values: the best
;; result, a boolean indicating whether or not it converged to within
;; tolerance, and the number of iterations that were attempted. These
;; macros bind _ to the numerical result and @ to the number of
;; iterations. They work like their "normal" namesakes if used on
;; other (i.e. non-library) functions.

(defmacro if_ (if then &optional else)
  "make a decision based on whether or not `IF` converged and execute `THEN` or `ELSE`
accordingly, with `_` bound to the result and `@` bound to the number of attempts."
  (let ((l (gensym)))
    `(let ((,l (multiple-value-call #'list ,if)))
       (if (eq 3 (length ,l))
           (let ((_ (car ,l))
                 (@ (caddr ,l)))
             (declare (ignorable @ _))
             (if (cadr ,l)
                 ,then
                 ,else))
           (let ((_ (car ,l))
                 (@ nil))
             (declare (ignorable @))
             (if _
                 ,then
                 ,else))))))

(defmacro when_ (if then)
  "Execute `THEN` if `IF` converges, with `_` bound to the result and `@`
bound to the number of attempts."
  `(if_ ,if ,then nil))

(defmacro cond_ (&rest clauses)
  "Execute clauses based on whether or not functions converge, with `_` bound
to the result and `@` bound to the number of attempts."
  (when clauses
    (let ((cl1 (car clauses)))
      `(if_ ,(car cl1)
            ,@(cdr cl1)
            (cond_ ,@(cdr clauses))))))

;; "safe" version of divide to use internally
(defun /_ (a &rest b)
  "divide, but avoid singularities"
  (let ((b (reduce #'* b)))
    (/ a (if (zerop b) +epsilon+ b))))

;; define some utilities for stuffing expressions into lambdas as necessary
  
(defun free-vars (e)
  (labels ((rec (e)
             (if (atom e)
                 (when (and
                        (not (null e))
                        (symbolp e)
                        (not (boundp e))
                        (not (fboundp e)))
                   (list e))
                 (append (rec (car e)) (rec (cdr e))))))
    (remove-duplicates (rec e) :from-end t)))
  
(defmacro genlambda/1 (e)
  "take an expression and generate a 1-argument lambda"
  (if (and (consp e)
           (or
            (eq (car e) 'FUNCTION)
            (eq (car e) 'LAMBDA)))
      `,e
      (let ((fv (free-vars e)))
        (cond ((eq 1 (length fv)) `(lambda ,fv ,e))
              ((eq 0 (length fv)) (if (atom e)
                                      (if (null e)
                                          `(lambda (x) (declare (ignore x)) nil)
                                          `(lambda (x) (,e x)))
                                      `(lambda (r) (apply (quote ,(car e)) (append (list ,@(cdr e)) (list r))))))
              (t (error "too many free variables - should have 1 at most"))))))

(defmacro genlambda/2 (e)
  "take an expression and generate a 2-argument lambda"
  (if (and (consp e)
           (or
            (eq (car e) 'FUNCTION)
            (eq (car e) 'LAMBDA)))
      `,e
      (let ((fv (free-vars e)))
        (cond ((eq 2 (length fv)) `(lambda ,fv ,e))
              ((eq 1 (length fv))
               (if (atom e)
                   `(lambda (x y) (,e x y))
                   (let ((x (car fv))
                         (y (gensym)))
                     `(lambda (,x ,y) (declare (ignore ,y)) ,e))))
              ((eq 0 (length fv))
               `(lambda (&rest r) (apply (quote ,(car e)) (append (list ,@(cdr e)) r))))
              (t (error "too many free variables - should have 2 at most"))))))

