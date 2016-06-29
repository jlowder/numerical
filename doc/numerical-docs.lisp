(defpackage :numerical-docs
  (:use #:cl #:gendoc
        :numerical
        :numerical.integration 
        :numerical.roots
        :numerical.diffeq
        :numerical.interpolation)
  (:export #:generate))

(in-package :numerical-docs)

(defun generate ()
  (gendoc (:output-filename "ref.html"
                            :css "simple.css")
          (:apiref :numerical
                   :numerical.integration 
                   :numerical.roots
                   :numerical.interpolation)))
