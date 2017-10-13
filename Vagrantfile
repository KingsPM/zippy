# -*- mode: ruby -*-
# vi: set ft=ruby :
Vagrant.configure(2) do |config|
  # share type dependent on platform
  host = RbConfig::CONFIG['host_os']
  sharetype = host =~ /darwin/ || host =~ /linux/ ? "nfs" : "smb"
  config.vm.synced_folder ".", "/home/vagrant/zippy", type: sharetype
  config.vm.box = "ubuntu/precise64"
  # port config
  config.vm.network "forwarded_port", guest: 80, host: 5000
  config.vm.network "forwarded_port", guest: 5000, host: 5001
  config.vm.network "private_network", ip: "55.55.55.5"
  # proxy config (CNTLM)
  if Vagrant.has_plugin?("vagrant-proxyconf")
    config.env_proxy.http  = "http://10.0.2.2:3128"
    config.env_proxy.https = "http://10.0.2.2:3128"
    config.proxy.http      = "http://10.0.2.2:3128"
    config.proxy.https     = "http://10.0.2.2:3128"
    config.proxy.no_proxy  = "localhost,127.0.0.1,192.168.33.*"
  end
  config.vm.provider "virtualbox" do |vb|
    vb.gui = false
    vb.cpus = 2
    vb.memory = 6144
  end
end
