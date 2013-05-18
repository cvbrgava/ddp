#!/usr/bin/expect -f

# connect via scp
spawn-fcgi scp html.tar.gz ee08b037@10.7.0.3:~/
#######################
expect {
-re ".*es.*o.*" {
exp_send "yes\r"
exp_continue
}
-re ".*sword.*" {
exp_send "cvbrgava\r"
}
}
interact
