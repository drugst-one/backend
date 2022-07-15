# Generated by Django 3.0.5 on 2020-04-03 16:36

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('drugstone', '0025_auto_20200403_1457'),
    ]

    operations = [
        migrations.AddField(
            model_name='drug',
            name='drug_id',
            field=models.CharField(default='No id', max_length=10, unique=True),
        ),
        migrations.AlterField(
            model_name='drug',
            name='name',
            field=models.CharField(default='No name', max_length=128),
        ),
        migrations.AlterField(
            model_name='drug',
            name='status',
            field=models.CharField(default='No status', max_length=32),
        ),
    ]