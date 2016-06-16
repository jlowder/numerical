(defsystem :numerical-docs
  :version "0.1.0"
  :author "Jason Lowdermilk <jlowdermilk@gmail.com>"
  :description "Document generation for numerical"
  :depends-on (:cl-gendoc :numerical)
  :license "MIT"
  :serial t
  :components
  ((:file "numerical-docs")))

(defmethod perform :after ((o load-op) (c (eql (find-system :numerical-docs))))
  (let ((fn (find-symbol (symbol-name 'generate) (find-package :numerical-docs))))
    (funcall fn)))

(defmethod operation-done-p ((o load-op) (c (eql (find-system :numerical-docs))))
  nil)
