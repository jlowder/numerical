(defsystem numerical-test
  :author "Jason Lowdermilk <jlowdermilk@gmail.com>"
  :license "MIT"
  :description "Unit tests for numerical"
  :depends-on (:lisp-unit :numerical)
  :components ((:file "test")))
