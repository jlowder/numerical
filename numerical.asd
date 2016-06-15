(defsystem numerical
    :name "numerical"
    :version "0.1.0"
    :author "Jason Lowdermilk <jlowdermilk@gmail.com>"
    :license "MIT"
    :description "Library of numerical methods"
    :serial t
    :components
    ((:module "src"
              :components
              ((:file "numerical")
               (:file "roots")
;;               (:file "interpolation")
;;               (:file "integration")))))
                        
))))
