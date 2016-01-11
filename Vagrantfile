# -*- mode: ruby -*-
# vi: set ft=ruby :
Vagrant.configure(2) do |config|
  config.vm.synced_folder ".", "/home/vagrant/dev/zippy" #, type: "nfs"
  config.vm.box = "ubuntu/precise64"
  config.vm.network "forwarded_port", guest: 8000, host: 8000
  config.vm.network "forwarded_port", guest: 8080, host: 8080
  config.vm.network "forwarded_port", guest: 5000, host: 5000
  config.vm.network "private_network", ip: "55.55.55.5"
  config.vm.provider "virtualbox" do |vb|
    vb.gui = false
    vb.cpus = 2
    vb.memory = 8192
  end
end
